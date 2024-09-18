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

# 定义肽段的单电荷状态的单同位素质量（从PeptideMass获取）
monoisotopic_mass <- 785.3247

# 定义要检查的电荷状态范围
charge_states <- 1:5  # 检查所有电荷状态

# 计算不同电荷状态下的 m/z 值
mz_values <- sapply(charge_states, function(z) {
  (monoisotopic_mass + (z * 1.007276)) / z
})

# 输出计算的 m/z 值，调试信息
cat("Calculated m/z values:\n")
print(mz_values)

# 设置 m/z 提取的容差范围
mz_tolerance <- 0.5  # 尝试增加容差范围

# 初始化一个列表来存储每个电荷状态的 EIC 结果
eic_results <- list()

# 遍历每个 m/z 值并逐个提取离子色谱图（EICs）
for (i in seq_along(mz_values)) {
  mz <- mz_values[i]
  cat("Extracting chromatogram for m/z:", mz, "with tolerance:", mz_tolerance, "\n")
  
  # 提取 EIC 并添加错误处理
  eic <- tryCatch({
    chromatogram(raw_data, mz = c(mz - mz_tolerance, mz + mz_tolerance))
  }, error = function(e) {
    message("Error extracting chromatogram for m/z ", mz, ": ", e)
    NULL  # 返回NULL以继续处理其他 m/z
  })
  
  # 检查提取的结果并添加到列表中
  if (!is.null(eic)) {
    eic_results[[i]] <- eic
  } else {
    eic_results[[i]] <- NULL
    cat("No data found for m/z:", mz, "\n")
  }
}

# 过滤掉失败的提取
eic_results <- eic_results[!sapply(eic_results, is.null)]

# 输出实际提取到的 EIC数
cat("Number of EICs extracted:", length(eic_results), "\n")

# 自动选择有信号的 EIC
eic_with_signal <- list()

for (eic in eic_results) {
  # 使用 lapply 方法检查每个色谱图中的信号
  has_signal <- lapply(eic, function(x) {
    any(intensity(x) > 0, na.rm = TRUE)  # 忽略缺失值
  })
  
  if (any(unlist(has_signal), na.rm = TRUE)) {  # 如果有信号，忽略缺失值
    eic_with_signal[[length(eic_with_signal) + 1]] <- eic
  }
}

# 检查是否找到有信号的 EIC
if (length(eic_with_signal) > 0) {
  cat("Found", length(eic_with_signal), "EIC(s) with signal.\n")
} else {
  cat("No EICs with signal were found.\n")
}


# Output detailed information of EICs with signal
for (i in seq_along(eic_with_signal)) {
  eic <- eic_with_signal[[i]]
  
  cat("EIC", i, ":\n")
  
  # Iterate through each chromatogram in the EIC
  for (chrom in eic) {
    # Extract m/z range and retention time information
    mz_range <- range(chrom@mz)
    rt_range <- range(rtime(chrom))
    
    # Print the extracted information
    cat("  m/z range:", mz_range, "\n")
    cat("  Retention time range:", rt_range, "\n")
    cat("  Number of points:", length(rtime(chrom)), "\n")
  }
}

# 提取 m/z = 786.332 和 393.6696 的 EIC 数据
target_mz_786 <- 786.332
target_mz_394 <- 393.6696

# 检查目标 m/z 是否存在于计算的 m/z 值中
index_786 <- which(abs(mz_values - target_mz_786) < 0.001)
index_394 <- which(abs(mz_values - target_mz_394) < 0.001)

# 检查并可视化EIC for m/z 786.332
if (length(index_786) > 0) {
  eic_786 <- eic_results[[index_786]]
  if (!is.null(eic_786)) {
    eic_786_df <- data.frame(
      rt = eic_786[[1]]@rtime,
      intensity = eic_786[[1]]@intensity
    )
    
    ggplot(eic_786_df, aes(x = rt, y = intensity)) +
      geom_line(color = "blue") +
      theme_minimal() +
      labs(title = "EIC for m/z: 786.332",
           x = "Retention Time (seconds)",
           y = "Intensity") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red")
  } else {
    cat("No EIC found for m/z: 786.332\n")
  }
} else {
  cat("No EIC found for m/z: 786.332\n")
}

# 检查并可视化EIC for m/z 393.6696
if (length(index_394) > 0) {
  eic_394 <- eic_results[[index_394]]
  if (!is.null(eic_394)) {
    eic_394_df <- data.frame(
      rt = eic_394[[1]]@rtime,
      intensity = eic_394[[1]]@intensity
    )
    
    ggplot(eic_394_df, aes(x = rt, y = intensity)) +
      geom_line(color = "Blue") +
      theme_minimal() +
      labs(title = "EIC for m/z: 393.6696",
           x = "Retention Time (seconds)",
           y = "Intensity") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red")
  } else {
    cat("No EIC found for m/z: 393.6696\n")
  }
} else {
  cat("No EIC found for m/z: 393.6696\n")
}



