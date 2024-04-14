library(Seurat)
add_score <- function(seurat_obj, gene_list) {
  # 检查Seurat对象和基因列表是否有效
  if (is.null(seurat_obj) || !is(seurat_obj, "Seurat")) {
    stop("Invalid Seurat object provided.")
  }
  
  if (is.null(gene_list) || !is.list(gene_list)) {
    stop("Invalid gene list provided.")
  }
  
  # 对于基因列表中的每个细胞类型
  for (cell_type in names(gene_list)) {
    genes <- gene_list[[cell_type]]
    
    # 检查基因是否在Seurat对象中
    if (!all(genes %in% rownames(seurat_obj))) {
      warning(paste("Some genes for", cell_type, "are not found in the Seurat object. Skipping this cell type."))
      
    }
    
    # 计算并添加模块分数
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(genes),name=cell_type
      
    )
  }
  return(seurat_obj)
}
