# SeuratQCUtils
# Takes Seurat object
# Package 1 ---------------------------------------------------------------
addQCMetrics_seu <- function(seu){
  ## Per cell
  seu$percent.mt <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  libsize_drop <- isOutlier(
    metric = as.numeric(unlist(seu@meta.data[, grepl("nCount", colnames(seu@meta.data))])), 
    type = "lower",
    log = TRUE) 
  
  mito_drop <- isOutlier(
    metric = seu$percent.mt,
    type = "higher") 
  
  seu$libsize_drop <- libsize_drop
  seu$mito_drop <- mito_drop
  
  ## Per gene
  gene_means <- as.numeric(unlist(rowMeans(GetAssayData(seu, "counts"), na.rm = TRUE)))
  lowgenecount_drop <- log(gene_means) < -5 | gene_means <= 0
  
  seu[[names(seu@assays)[1]]][["means"]] <- gene_means
  seu[[names(seu@assays)[1]]][["lowgenecount_drop"]] <- lowgenecount_drop
  
  return(seu)
}

addQCMetricsPerCell_seu <- function(seu){
  ## Per cell
  seu$percent.mt <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  libsize_drop <- isOutlier(
    metric = as.numeric(unlist(seu@meta.data[, grepl("nCount", colnames(seu@meta.data))])), 
    type = "lower",
    log = TRUE) 
  
  mito_drop <- isOutlier(
    metric = seu$percent.mt,
    type = "higher") 
  
  seu$libsize_drop <- libsize_drop
  seu$mito_drop <- mito_drop
  
  return(seu)
}

addQCMetricsPerGene_seu <- function(seu){
  ## Per gene
  gene_means <- as.numeric(unlist(rowMeans(GetAssayData(seu, "counts"), na.rm = TRUE)))
  lowgenecount_drop <- log(gene_means) < -5 | gene_means <= 0
  
  seu[[names(seu@assays)[1]]][["means"]] <- gene_means
  seu[[names(seu@assays)[1]]][["lowgenecount_drop"]] <- lowgenecount_drop
  
  return(seu)
}


plotHist <- function(plot_df, xvar = "logtotal", yvar = "libsize_drop"){
  p <- plot_df %>%
    ggplot(aes(x = get(xvar), fill = get(yvar))) +
    geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity', bins = 30) +
    scale_fill_manual(values = c("#404080", "#69b3a2")) +
    theme_ipsum() +
    labs(fill = "") +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 12))
  
  p
}

plot_Hist_Low_Lib_Sizes <- function(seu, yvar = "libsize_drop"){
  CD <- seu@meta.data
  CD[[yvar]] <- factor(CD[[yvar]])
  CD[[yvar]] <- relevel(CD[[yvar]], "TRUE")
  
  plot_df <- data.frame(logtotal = log(as.numeric(seu@meta.data[, grepl("nCount", colnames(seu@meta.data))])),
                        CD)
  
  p <- plotHist(plot_df, xvar = "logtotal", yvar = yvar)
  
  p <- p + xlab("Log cell level total count") + ylab("Frequency") + 
    ggtitle('Cells with low library size')
  
  p
}

plot_Hist_High_Mito_Props <- function(seu, yvar = "mito_drop"){
  CD <- seu@meta.data
  CD[[yvar]] <- factor(CD[[yvar]])
  CD[[yvar]] <- relevel(CD[[yvar]], "TRUE")
  
  plot_df <- data.frame(log_mito_percent = log(seu$percent.mt),
                        CD)
  
  p <- plotHist(plot_df, xvar = "log_mito_percent", yvar = yvar)
    
  p <- p + xlab("Log cell level mitochondria percent") + ylab("Frequency") + 
    ggtitle('Cells with high mitochondria percentage')
  
  p
}

plot_Hist_Low_Abun_Genes <- function(seu){
  plot_df <- data.frame(mean_genecount = log(unlist(seu[[names(seu@assays)[1]]][["means"]])),
                        lowgenecount_drop = factor(unlist(seu[[names(seu@assays)[1]]][["lowgenecount_drop"]])))
  
  plot_df$lowgenecount_drop <- relevel(plot_df$lowgenecount_drop, "TRUE")
  
  p <- plotHist(plot_df, xvar = "mean_genecount", yvar = "lowgenecount_drop")
  
  p <- p + xlab("Log mean count across all cells") + ylab("Frequency") + 
    ggtitle('Low abundance genes')
  
  p
}


# Spatial plot for binary flags 
plotQCSpatial_seu <- function(xe, flag = "libsize_drop"){
  if(!("array_row" %in% colnames(xe@meta.data))){
    xe[[c("array_row", "array_col")]] <- GetTissueCoordinates(xe)[, 1:2]
    xe[["cell_id"]] <- colnames(xe)
  }
  
  df_cellmeta <- xe@meta.data[, grepl("nCount_|nFeature_", colnames(xe@meta.data))]
  df <- data.frame(xe[[c("cell_id", "array_col", "array_row", flag)]])
  df <- cbind(df, df_cellmeta)
  ggplot(df, aes(x = array_row, y = array_col, color = get(flag))) + 
    geom_point(size = 0.3) + 
    coord_fixed() + 
    scale_color_manual(name = flag, values = c("gray85", "red")) + 
    ggtitle("QC dots") + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
}


plotBar_seu <- function(seu, decreasing = TRUE, ylabel = "Percentage", 
                             title = "Cell-type Frequency in Annotated Chromium"){
  CD <- seu@meta.data
  cnt <- plyr::count(CD$Annotation)
  CD$Annotation <- factor(CD$Annotation, 
                          levels = cnt$x[order(cnt$freq, decreasing = decreasing)])
  
  p <- ggplot(data = CD, aes(x = Annotation)) + 
    geom_bar(aes(y = (..count..)/sum(..count..))) +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    ylab(ylabel) + ggtitle(title)
  
  p
}


# plotOnUMAP
# Takes sce object
# Package Number 2 --------------------------------------------------------
plotFeatureUMAP <- function(sce, feat_name = "ADIPOQ", facet = NULL,
                            dim_red_name = "UMAP.HARMONY"){
  mat <- as.data.frame(t(logcounts(sce)))
  meta <- as.data.frame(colData(sce))
  dim.red <- as.data.frame(reducedDim(sce, dim_red_name))
  CD <- cbind(mat, meta, dim.red)
  
  p <- ggplot(CD, aes(x = UMAP1, y = UMAP2, color = get(feat_name))) +
    geom_point(size = 0.3) + 
    labs(color = NULL) + 
    theme_bw() + 
    theme(legend.position="right", panel.border = element_blank()) +
    ggtitle(feat_name) 
  
  if(!is.null(facet)){
    p <- p + facet_wrap(~get(facet))
  } 
  
  p
}

plotClusterUMAP <- function(sce, cluster_name = "spatial.cluster", facet = NULL, 
                            dim_red_name = "UMAP.HARMONY"){
  meta <- as.data.frame(colData(sce))
  dim.red <- as.data.frame(reducedDim(sce, dim_red_name))
  CD <- cbind(meta, dim.red)
  
  p <- ggplot(CD, aes(x = UMAP1, y = UMAP2, color = as.factor(get(cluster_name)))) +
    geom_point(size = 0.3) + 
    labs(color = NULL) + 
    theme_bw() + 
    theme(legend.position="right", panel.border = element_blank()) +
    scale_color_manual(name = cluster_name, 
                       values = scales::hue_pal()(length(unique(sce[[cluster_name]])))) +
    ggtitle(cluster_name) +
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  if(!is.null(facet)){
    p <- p + facet_wrap(~get(facet))
  } 
  
  p
}


# # Spatial plot for continuous features
# SpatialFeaturePlot_cont <- function(xe, color_by = "ABCC11"){
#   CD <- xe@meta.data
#   feat_mat <- GetAssayData(object = xe, slot = "counts")
#   CD <- cbind(CD, t(feat_mat))
#   
#   ggplot(CD,
#          aes(x = array_row, y = array_col,
#              color = get(color_by))) +
#     geom_point(size = 0.3) +
#     theme_void() + 
#     theme(legend.position="right") +
#     scale_colour_gradientn(name = color_by, colors = myPalette(100))
# }

# # Spatial plot for categorical value with self-defined palette
# SpatialFeaturePlot_cate <- function(xe, color_by = "seurat_clusters", palette = NULL){
#   CD <- xe@meta.data
#   
#   ggplot(CD,
#          aes(x = array_row, y = array_col,
#              color = get(color_by))) +
#     geom_point(size = 0.3) +
#     theme_void() + 
#     theme(legend.position="right") +
#     scale_color_manual(name = color_by, values = palette) +
#     guides(colour = guide_legend(override.aes = list(size = 3)))
# }
# 

# extractCol <- function(pumap){
#   group <- as.numeric(ggplot_build(pumap)$data[[1]]$group)
#   col <- ggplot_build(pumap)$data[[1]]$colour
#   col <- forcats::fct_reorder(col, group)
#   cols_Red <- levels(col)
#   
#   return(cols_Red)
# }
