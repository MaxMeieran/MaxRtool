clusterCornerAxes_my=function (object = NULL, reduction = "umap", groupFacet = groupFacet, 
          clusterCol = "seurat_clusters", customColors = NULL, pSize = 1, aspect.ratio = NULL, 
          noSplit = TRUE, nrow = 1, relLength = 0.25, relDist = 0.1, 
          axes = "mul", show.legend = TRUE, legendPos = "right", keySize = 5, 
          cellLabel = FALSE, cellLabelSize = 6, cellLabelColor = "black", 
          lineTextcol = "black", stripCol = "white", arrowType = "closed", 
          cornerTextSize = 3, base_size = 14, themebg = "default", 
          addCircle = FALSE, cicAlpha = 0.1, cicLineSize = 1, cicLineColor = "grey50", 
          cicLineLty = "dashed", nbin = 100, nsm = 10, addsm = 1, 
          qval = 1, sfac = 1.5) 
{
  reduc <- data.frame(Seurat::Embeddings(object, reduction = reduction))
  meta <- object@meta.data
  pc12 <- cbind(reduc, meta)
  namePos <- pc12 %>% dplyr::group_by(.data[[clusterCol]]) %>% 
    dplyr::summarise(posMedia1 = stats::median(get(colnames(pc12)[1])), 
                     posMedia2 = stats::median(get(colnames(pc12)[2])))
  range <- floor(min(min(pc12[, 1]), min(pc12[, 2])))
  lower <- range - relDist * abs(range)
  labelRel <- relDist * abs(lower)
  linelen <- abs(relLength * lower) + lower
  mid <- abs(relLength * lower)/2 + lower
  if (startsWith(reduction, "umap")) {
    axs_label <- paste("UMAP", 2:1, sep = "")
  }
  else if (startsWith(reduction, "tsne")) {
    axs_label <- paste("t-SNE", 2:1, sep = "")
  }
  else {
    print("Please give correct type(umap or tsne)!")
  }
  if (axes == "mul") {
    axes <- data.frame(x1 = c(lower, lower, lower, linelen), 
                       y1 = c(lower, linelen, lower, lower), linegrou = c(1, 
                                                                          1, 2, 2))
    label <- data.frame(lab = c(axs_label), angle = c(90, 
                                                      0), x1 = c(lower - labelRel, mid), y1 = c(mid, lower - 
                                                                                                  labelRel))
  }
  else if (axes == "one") {
    firstFacet <- unique(pc12[, groupFacet])[1]
    axes <- data.frame(x1 = c(lower, lower, lower, linelen), 
                       y1 = c(lower, linelen, lower, lower), linegrou = c(1, 
                                                                          1, 2, 2), group = rep(firstFacet, 2))
    label <- data.frame(lab = c(axs_label), angle = c(90, 
                                                      0), x1 = c(lower - labelRel, mid), y1 = c(mid, lower - 
                                                                                                  labelRel), group = rep(firstFacet, 2))
    colnames(axes)[4] <- groupFacet
    colnames(label)[5] <- groupFacet
  }
  else {
    print("Please give correct args(mul or one)!")
  }
  p <- ggplot2::ggplot(pc12, ggplot2::aes_string(x = colnames(pc12)[1], 
                                                 y = colnames(pc12)[2])) + ggplot2::geom_point(ggplot2::aes_string(color = clusterCol), 
                                                                                               size = pSize, show.legend = show.legend) + ggplot2::theme_classic(base_size = base_size) + 
    ggplot2::labs(x = "", y = "") + ggplot2::theme(strip.background = ggplot2::element_rect(colour = NA, 
                                                                                            fill = stripCol), aspect.ratio = aspect.ratio, legend.position = legendPos, 
                                                   plot.title = ggplot2::element_text(hjust = 0.5), axis.line = ggplot2::element_blank(), 
                                                   axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank()) + 
    ggplot2::geom_line(data = axes, ggplot2::aes(x = x1, 
                                                 y = y1, group = linegrou), color = lineTextcol, 
                       arrow = ggplot2::arrow(length = ggplot2::unit(0.1, 
                                                                     "inches"), ends = "last", type = arrowType)) + 
    ggplot2::geom_text(data = label, ggplot2::aes(x = x1, 
                                                  y = y1, angle = angle, label = lab), color = lineTextcol, 
                       fontface = "italic", size = cornerTextSize) + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = keySize)))
  if (!is.null(customColors)) {
    p <- p + ggplot2::scale_color_manual(values = customColors)
  }
  if (cellLabel == FALSE) {
    plabel <- p
  }
  else {
    plabel <- p + ggrepel::geom_text_repel(data = namePos, 
                                           ggplot2::aes_string(x = "posMedia1", y = "posMedia2", 
                                                               label = clusterCol), show.legend = F, size = cellLabelSize, 
                                           color = cellLabelColor)
  }
  if (addCircle == FALSE) {
    p0 <- plabel
  }
  else {
    p0 <- plabel + ggunchull::stat_unchull(ggplot2::aes_string(fill = clusterCol), 
                                           alpha = cicAlpha, size = cicLineSize, color = cicLineColor, 
                                           lty = cicLineLty, show.legend = F, nbin = nbin, 
                                           nsm = nsm, addsm = addsm, sfac = sfac, qval = qval)
  }
  if (noSplit == TRUE) {
    p1 <- p0
  }
  else {
    p1 <- p0 + ggplot2::facet_wrap(facets = groupFacet, 
                                   nrow = nrow)
  }
  if (themebg == "bwCorner") {
    p2 <- p1 + ggplot2::theme_bw(base_size = base_size) + 
      ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                     axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), 
                     aspect.ratio = 1, strip.background = ggplot2::element_rect(colour = NA, 
                                                                                fill = stripCol))
  }
  else if (themebg == "default") {
    p2 <- p1
  }
  # 应用自定义颜色（如果提供了）

  return(p2)
}
