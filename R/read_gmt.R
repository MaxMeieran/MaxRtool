read_gmt <- function(file_path) {
  # 打开文件
  con <- file(file_path, "r")
  
  # 初始化列表来存储基因集
  gene_sets <- list()
  
  # 逐行读取文件
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {  # 如果到达文件末尾，跳出循环
      break
    }
    # 分割每行数据
    elements <- strsplit(line, "\t")[[1]]
    set_name <- elements[1]  # 基因集的名称
    # 可能需要调整下标以适应是否有描述字段
    genes <- elements[-c(1,2)]  # 基因列表，假设第二项是描述
    gene_sets[[set_name]] <- genes
  }
  
  # 关闭文件连接
  close(con)
  
  # 返回基因集列表
  return(gene_sets)
}