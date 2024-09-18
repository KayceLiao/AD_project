# 加载必要的R包
library(xcms)
library(MSnbase)
library(ggplot2)

# 确保文件路径正确并加载文件
file_path <- "C:/Users/77977/OneDrive/Desktop/dataset/UWA558_B.mzXML"

# 检查文件是否存在
if (!file.exists(file_path)) {
  stop("The file does not exist at the specified path. Please check the path.")
}

# 加载 mzXML 文件
raw_data <- readMSData(files = file_path, mode = "onDisk")

# 使用 centWave 方法进行峰检测
cwp <- CentWaveParam(ppm = 15, peakwidth = c(5, 20), snthresh = 10)
xdata <- findChromPeaks(raw_data, param = cwp)

# 定义不同电荷状态的质量差
mass_differences <- 26.015650 / (1:5)  # 电荷状态 z = 1 到 5

# 定义时间窗口（秒）
time_window <- 5  # 可以根据需要调整

# 在质谱数据的一小部分中工作（保留时间窗口）
small_rt_window <- c(1000, 1500)  # 定义保留时间窗口
small_data <- filterRt(xdata, rt = small_rt_window)

# 提取选定时间窗口内的峰
peaks <- chromPeaks(small_data)

# 初始化一个空列表以存储匹配对
matching_pairs <- list()

# 遍历每对峰以找到匹配
for (i in 1:(nrow(peaks) - 1)) {
  for (j in (i + 1):nrow(peaks)) {
    # 计算时间差
    time_diff <- abs(peaks[i, "rt"] - peaks[j, "rt"])
    
    # 检查时间差是否在时间窗口内
    if (time_diff <= time_window) {
      # 计算 m/z 差异
      mz_diff <- abs(peaks[i, "mz"] - peaks[j, "mz"])
      
      # 检查 mz_diff 是否匹配任何质量差异
      for (z in 1:5) {
        if (abs(mz_diff - mass_differences[z]) < 0.01) {  # 容差可以根据需要调整
          matching_pairs[[length(matching_pairs) + 1]] <- list(
            peak1 = peaks[i, ],
            peak2 = peaks[j, ],
            charge = z
          )
        }
      }
    }
  }
}

# 检查是否找到匹配对
if (length(matching_pairs) > 0) {
  # 打印找到的匹配对
  cat("Found matching pairs:\n")
  for (pair in matching_pairs) {
    cat("Peak 1 m/z:", pair$peak1["mz"], "RT:", pair$peak1["rt"],
        "Peak 2 m/z:", pair$peak2["mz"], "RT:", pair$peak2["rt"],
        "Charge:", pair$charge, "\n")
  }
  
  # 提取匹配对信息用于可视化
  match_data <- data.frame(
    peak1_mz = sapply(matching_pairs, function(x) x$peak1["mz"]),
    peak1_rt = sapply(matching_pairs, function(x) x$peak1["rt"]),
    peak2_mz = sapply(matching_pairs, function(x) x$peak2["mz"]),
    peak2_rt = sapply(matching_pairs, function(x) x$peak2["rt"]),
    charge = sapply(matching_pairs, function(x) x$charge)
  )
  
  # 可视化峰对
  p <- ggplot(match_data) +
    geom_point(aes(x = peak1_rt, y = peak1_mz, color = as.factor(charge)), shape = 16, size = 2) +
    geom_point(aes(x = peak2_rt, y = peak2_mz, color = as.factor(charge)), shape = 16, size = 2) +
    geom_segment(aes(x = peak1_rt, y = peak1_mz, xend = peak2_rt, yend = peak2_mz, color = as.factor(charge)), linetype = "dashed") +
    scale_color_manual(values = c("red", "green", "blue", "orange", "purple"), name = "Charge State") +
    theme_minimal() +
    labs(title = "Detected Peak Pairs with Matching m/z Differences",
         x = "Retention Time (seconds)",
         y = "m/z")
  
  print(p)
} else {
  cat("No matching pairs found in the specified window.\n")
}
