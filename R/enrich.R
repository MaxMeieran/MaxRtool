library(clusterProfiler)
library(stringr)
library(org.Hs.eg.db)
library(dplyr)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(forcats)
library(deeptime)
#BiocManager::install("dittoSeq")
library(dittoSeq)

enrich=function(module,gene,ont='BP'){
  genelist <- bitr(gene, 
                   fromType="SYMBOL",
                   toType="ENTREZID", 
                   OrgDb='org.Hs.eg.db')
  genelist <- pull(genelist, ENTREZID)
  
  if (module =='go'){
    df_GO<- enrichGO(
      genelist, 
      OrgDb="org.Hs.eg.db",
      keyType='ENTREZID',
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = T
    )
  df_go<- data.frame(df_GO)
  df_go$GeneRatio <- sapply(strsplit(df_go$GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  df_go <- df_go%>% mutate(Plog = -log10(pvalue))
  return(df_go)}
  
  if (module =='go'){
    df_GO<- enrichGO(
      genelist, 
      OrgDb="org.Hs.eg.db",
      keyType='ENTREZID',
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = T
    )
  df_go<- data.frame(df_GO)
  df_go$GeneRatio <- sapply(strsplit(df_go$GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  df_go <- df_go%>% mutate(Plog = -log10(pvalue))
  return(df_go)}
  
  
  if (module=='kegg'){
    KEGG<-enrichKEGG(genelist,#KEGG富集分析
                       organism = 'hsa',
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
    kegg=as.data.frame(KEGG)
    ek.rt=kegg
    ek.rt <- separate(data=ek.rt, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
    ek.rt <- separate(data=ek.rt, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
    ek.rt <- mutate(ek.rt, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor
    return(kegg)}
}

enrichwithplot <- function(module, gene) {ont = 'ALL'
  # 将基因符号转换为ENTREZ ID
  genelist <- bitr(gene,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  genelist <- pull(genelist, ENTREZID)
  
  if (module == 'go') {
    # 进行GO富集分析
    df_GO <- enrichGO(
      gene          = genelist,
      OrgDb         = org.Hs.eg.db,
      keyType       = 'ENTREZID',
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    
    # 处理得到的数据框
    df_go <- as.data.frame(df_GO)
    df_go$GeneRatio <- sapply(strsplit(df_go$GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    df_go <- df_go %>% mutate(Plog = -log10(pvalue))
    
    # 提取前十名的富集结果
    top_hits <- df_go %>%
      arrange(desc(Plog)) %>%
      group_by(ONTOLOGY) %>%
      slice_head(n = 10)
    
    ontology_counts <- table(top_hits$ONTOLOGY)
    
    # 为每个类别指定一个颜色
    colors <- c("#BB3DFB", "#F90107", "#F88604")  # 确保这个颜色列表的顺序与ontology_counts的标签顺序相匹配
    
    # 根据每个类别的数量重复对应的颜色
    col <- rep(colors, times = ontology_counts)
    
    top_hits <- top_hits %>%
      arrange(ONTOLOGY, desc(Plog))
    
    # 使 Description 成为一个因子，水平按 Plog 的排序
    top_hits$Description <- factor(
      top_hits$Description, 
      levels = unique(top_hits$Description[order(top_hits$ONTOLOGY, -top_hits$Plog)])
    )
    
    # 绘制图表
    pp <- ggplot(top_hits, aes(x = Description, y = Plog, fill = ONTOLOGY)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(x = '', y = '-log10(p-value)', title = 'GO Enrichment Analysis') +
      scale_fill_manual(values = colors) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = 'none',
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.text.y = element_text(size = 12, colour = col, face = "bold"),
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))
    
    return(list(df_go = df_go, plot = pp))
  }

  if (module=='kegg'){
    print('KEGG还不支持！！！！！')
    }
}


enrich_gsego=function(module,DEG,ont='BP'){
  DEG<- DEG %>%
    arrange(desc(avg_log2FC))
  
  # 提取avg_log2FC列为向量
  fc_vector <- DEG$avg_log2FC
  fc_vector <- setNames(DEG$avg_log2FC, DEG$gene)
  
    # 使用 bitr 转换基因名
    genelist <- bitr(names(fc_vector), 
                     fromType = "SYMBOL",
                     toType = "ENTREZID", 
                     OrgDb = 'org.Hs.eg.db')
    
    # 过滤出成功转换的基因
    successful_conversions <- genelist[!is.na(genelist$ENTREZID), ]
    
    # 提取 ENTREZID 列，并设置为 avg_log2FC 向量的新名字
    valid_entrez_ids <- setNames(fc_vector[names(fc_vector) %in% successful_conversions$SYMBOL], 
                                 successful_conversions$ENTREZID)
    
    # 打印或返回 valid_entrez_ids
    print(valid_entrez_ids)
  
  
    df_GO<-  gseGO(
      valid_entrez_ids, 
      OrgDb="org.Hs.eg.db",
      keyType='ENTREZID',
      ont = 'ALL',
      minGSSize    = 100,
      maxGSSize    = 500,pvalueCutoff = 1,
      verbose      = FALSE
    
    )
  
  df_GO@result <- df_GO@result%>% mutate(Plog = -log10(p.adjust))
  return(df_GO)}
  
shuangce_go=function(frame,title='title'){
  if (nrow(frame) > 20) {
    top_max <- frame %>%
      filter(GeneRatio > 0)%>%
      arrange(desc(abs(GeneRatio))) %>%
      slice_head(n = 10)
    
    # 选择 GeneRatio 最小的10个且为负数
    top_min <- frame %>%
      filter(GeneRatio < 0) %>%
      arrange(GeneRatio) %>%
      slice_head(n = 10)
    
    # 合并结果
    frame <- rbind(top_max, top_min)
  }
  color <- rep("#ae4531", nrow(frame))
  color[which(frame$GeneRatio < 0)] <- "#2f73bb"
  
  frame$color <- color
  
  
  # 最终绘图代码：
ggplot(frame)+
    geom_col(aes(reorder(Description, GeneRatio), GeneRatio, fill = color))+
    scale_fill_manual(values = c("#2f73bb","#ae4531"))+
    # 加一条竖线：
    geom_segment(aes(y = 0, yend = 0,x = 0, xend = 21))+
    theme_classic()+
    ylim(-0.5,0.5)+
    coord_flip()+
    # 调整主题：
    theme(
      # 去除图例：
      legend.position = "none",
      # 标题居中：
      plot.title = element_text(hjust = 0.5),
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
    )+
    ylab("GeneRatio")+
    # 添加label：
  geom_text(data=data[which(data$value > 0), ],aes(x = pathway, y = 0, label = pathway), 
            hjust = 1.1, size = 4)+
  geom_text(data=data[which(data$value < 0), ],aes(x = pathway, y = 0, label = pathway), 
            hjust = -0.1, size = 4)+
    ggtitle(title)+
    scale_x_discrete(expand=expansion(add=c(0,1.5)))+
    # 添加箭头注释：
  geom_segment(aes(y = -0.03, yend = -0.35,x = 21, xend = 21),
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), 
               size = 0.5)+
  geom_segment(aes(y = 0.03, yend = 0.35,x = 21, xend = 21),
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), 
               size = 0.5)+
    annotate("text", x = 21, y = -0.5, label = "DOWN")+
    annotate("text", x = 21, y = 0.5, label = "UP")
  
}
   