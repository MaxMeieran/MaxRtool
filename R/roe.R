library(dplyr)

roe <- function(meta.data, group, celltype) {###group分组列名   celltype细胞注释列名
  # 全体样本中每种细胞类型的比例
  overall_prop <- prop.table(table(meta.data[[celltype]]))
  
  # 初始化一个数据框来存储结果
  results <- data.frame(Group = character(),
                        CellType = character(),
                        ROE = numeric(),
                        stringsAsFactors = FALSE)
  
  # 遍历每个分组
  unique_groups <- unique(meta.data[[group]])
  for (g in unique_groups) {
    # 提取该分组的数据
    group_data <- meta.data[meta.data[[group]] == g, ]
    
    # 计算该分组中每种细胞类型的比例
    group_prop <- prop.table(table(group_data[[celltype]]))
    
    # 计算ROE：分组比例 / 全体比例
    roe_values <- group_prop / overall_prop
    
    # 准备添加到结果数据框
    group_results <- data.frame(Group = rep(g, length(roe_values)),
                                CellType = names(roe_values),
                                ROE = as.numeric(roe_values),
                                stringsAsFactors = FALSE)
    
    # 将本组结果添加到总结果中
    results <- rbind(results, group_results)
  }
  
  # 返回结果数据框
  return(results)
}



