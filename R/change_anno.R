change_anno <- function(object, barcode, annocol, towhat) {
  # 检查barcode是否全部匹配
  column_names_in_row_names <- barcode %in% rownames(object@meta.data)
  all_true <- all(column_names_in_row_names)
  
  # 打印结果
  print(all_true)
  if (!all_true) {
    unmatched_barcodes <- barcode[!column_names_in_row_names]
    print("Unmatched barcodes:")
    print(unmatched_barcodes)
  } else {
    print("All barcodes are matched.")
  }
  
  # 仅在所有条码匹配时才执行更新
  if (all_true) {
    # 使用rlang::sym转换字符串为符号，用于动态列名引用
    annocol_sym <- rlang::sym(annocol)
    
    # 修改特定的列
    object@meta.data <- object@meta.data %>%
      mutate(!!annocol_sym := ifelse(row.names(.) %in% barcode, towhat, !!annocol_sym))
  }
  
  return(object)
}
