# 加载必要的R包
library(xcms)
library(MSnbase)
library(ggplot2)
library(plotly)
library(pracma)  

# 1. 读取数据并优化峰检测参数

# 加载 mzXML 文件
file_path <- "C:/Users/77977/OneDrive/Desktop/dataset/UWA558_B.mzXML"
raw_data <- readMSData(files = file_path, mode = "onDisk")


# 定义最佳参数组合
best_params <- list(ppm = 25, peakwidth = c(10, 50), snthresh = 5)

# 使用最佳参数组合进行峰检测
cwp <- CentWaveParam(ppm = best_params$ppm, peakwidth = best_params$peakwidth, snthresh = best_params$snthresh)
best_xdata <- findChromPeaks(raw_data, param = cwp)


# 2. 自相关分析

# 提取峰强度数据
chromatograms <- chromatogram(best_xdata)

# 初始化强度数据向量
intensity <- c()

# 检查并提取每个色谱图的强度
for (i in seq_along(chromatograms)) {
  ch <- chromatograms[[i]]  # 获取单个色谱图
  if (length(ch@intensity) > 0) {  # 确保有强度数据
    intensity <- c(intensity, ch@intensity)  # 合并强度数据
  }
}

# 打印强度数据的前几个值
cat("Intensity data (first few values):", head(intensity), "\n")

# 检查强度数据的长度
if (length(intensity) > 1) {
  # 计算自相关函数
  acf_values <- acf(intensity, plot = FALSE)$acf
  
  # 找到自相关峰（周期性）
  peak_positions <- which(acf_values > 0.5)  # 自定义阈值选择自相关峰
  cat("ACF Peak positions:", peak_positions, "\n")
} else {
  warning("强度数据长度不足，无法进行自相关分析")
}
acf(intensity, main = "ACF of Intensity Data")

# 3. 滑动窗口归一化

# 定义滑动窗口大小
window_size <- 50

# 使用滑动窗口计算平均值作为归一化因子
normalized_intensity <- sapply(1:(length(intensity) - window_size + 1), function(i) {
  mean(intensity[i:(i + window_size - 1)])  # 计算窗口内的平均值
})

# 确保 intensity 和 normalized_intensity 的长度匹配
# 使用向量化操作时，确保强度数据与归一化因子的长度一致
intensity_adjusted <- intensity[1:length(normalized_intensity)]  # 截取匹配的长度

# 使用归一化因子对调整后的强度数据进行归一化
intensity_normalized <- intensity_adjusted / normalized_intensity

# 打印归一化后的强度数据的前几个值
cat("Normalized Intensity data (first few values):", head(intensity_normalized), "\n")


# 创建数据框用于可视化
intensity_df <- data.frame(
  Index = 1:length(intensity_normalized),
  NormalizedIntensity = intensity_normalized
)

# 使用 ggplot2 绘制归一化后的强度数据
ggplot(intensity_df, aes(x = Index, y = NormalizedIntensity)) +
  geom_line(color = 'blue') +  # 绘制折线图
  theme_minimal() +
  labs(title = "Normalized Intensity Over Time",
       x = "Index",
       y = "Normalized Intensity") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")  # 在 y=1 处添加参考线

# 4. 峰对匹配算法

# 定义质量差和时间窗口
mass_differences <- list(`26.015650` = 26.015650 / (1:5), `28.031300` = 28.031300 / (1:5))
time_window <- 10
mass_tolerance <- 0.2

# 扩大保留时间窗口
small_rt_window <- c(900, 1600)
small_data <- filterRt(best_xdata, rt = small_rt_window)

# 提取峰数据
peaks <- chromPeaks(small_data)

# 检查 peaks 是否为空
if (nrow(peaks) == 0) {
  stop("No valid peaks found.")
}

# 转换峰数据为数据框
peaks_df <- as.data.frame(peaks)

# 初始化一个空数据框以存储匹配对
matching_pairs <- data.frame(
  peak1_mz = numeric(0), peak1_rt = numeric(0),
  peak2_mz = numeric(0), peak2_rt = numeric(0),
  charge = integer(0), mass_diff_type = character(0)
)

# 遍历找到匹配的峰对
for (i in 1:(nrow(peaks_df) - 1)) {
  for (j in (i + 1):nrow(peaks_df)) {
    time_diff <- abs(peaks_df[i, "rt"] - peaks_df[j, "rt"])
    
    if (is.na(time_diff) || time_diff > time_window) next
    
    mz_diff <- abs(peaks_df[i, "mz"] - peaks_df[j, "mz"])
    
    for (z in 1:5) {
      if (abs(mz_diff - mass_differences$`26.015650`[z]) < mass_tolerance) {
        matching_pairs <- rbind(matching_pairs, data.frame(
          peak1_mz = peaks_df[i, "mz"], peak1_rt = peaks_df[i, "rt"],
          peak2_mz = peaks_df[j, "mz"], peak2_rt = peaks_df[j, "rt"],
          charge = z, mass_diff_type = "26.015650"
        ))
      } else if (abs(mz_diff - mass_differences$`28.031300`[z]) < mass_tolerance) {
        matching_pairs <- rbind(matching_pairs, data.frame(
          peak1_mz = peaks_df[i, "mz"], peak1_rt = peaks_df[i, "rt"],
          peak2_mz = peaks_df[j, "mz"], peak2_rt = peaks_df[j, "rt"],
          charge = z, mass_diff_type = "28.031300"
        ))
      }
    }
  }
}

# 输出找到的匹配对
if (nrow(matching_pairs) > 0) {
  cat("Found matching pairs:\n")
  print(matching_pairs)
} else {
  cat("No matching pairs found in the specified window.\n")
}

# 5. 数据可视化

# 5.1 使用 ggplot2 绘制峰分布
ggplot(peaks_df, aes(x = rt, y = mz)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(x = "Retention Time (RT)", y = "m/z", title = "Detected Peaks Distribution")

# 5.2 使用 plotly 绘制3D可视化
plot_ly(peaks_df, x = ~rt, y = ~mz, z = ~into, type = 'scatter3d', mode = 'markers') %>%
  layout(scene = list(xaxis = list(title = 'Retention Time (RT)'),
                      yaxis = list(title = 'm/z'),
                      zaxis = list(title = 'Intensity')))
