df <- data.frame(
  gene_mean_vis = rowMeans(counts(vis_sce)), 
  gene_mean_aggxe = rowMeans(counts(aggxe_sce))
)
ggplot(df, aes(x = log(gene_mean_vis), y = log(gene_mean_aggxe))) +
  geom_point(size = 2) + 
  geom_abline()

# Check gene wise corr 

get_cor <- function(a_mat, b_mat,
                    xwise_cor = "cell-wise"){
  stopifnot(xwise_cor %in% c("cell-wise", "gene-wise"))
  
  ## Transpose
  if(xwise_cor == "gene-wise"){
    a_mat <- a_mat |> t()
    b_mat <- b_mat |> t()
  }
  
  cor_vector <- mapply(cor, 
                       as.data.frame(a_mat), 
                       as.data.frame(b_mat), 
                       use = "na.or.complete", 
                       method = "spearman")
  
  return(cor_vector)
}

vis_sce_mat <- counts(vis_sce)
aggxe_sce_mat <- counts(aggxe_sce)

# Order the spots and genes to be the same
aggxe_sce_mat <- aggxe_sce_mat[rownames(vis_sce_mat), colnames(vis_sce_mat)]



# gene_wise_cor_vis_aggxe <- get_cor(vis_sce_mat, aggxe_sce_mat, xwise_cor = "gene-wise")
gene_wise_cor_vis_aggxe <- get_cor(counts(vis_sce), counts(aggxe_sce), xwise_cor = "gene-wise")
plot(gene_wise_cor_vis_aggxe)
summary(gene_wise_cor_vis_aggxe)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.02485  0.12908  0.23270  0.27907  0.42898  0.69490 

# spot_wise_cor_vis_aggxe <- get_cor(vis_sce_mat, aggxe_sce_mat, xwise_cor = "cell-wise")
spot_wise_cor_vis_aggxe <- get_cor(counts(vis_sce), counts(aggxe_sce), xwise_cor = "cell-wise")
plot(spot_wise_cor_vis_aggxe)
summary(spot_wise_cor_vis_aggxe)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03327 0.42627 0.51596 0.52720 0.64088 0.81426 

df <- data.frame(cor = c(gene_wise_cor_vis_aggxe, spot_wise_cor_vis_aggxe), 
                 type = c(rep("gene_wise", length(gene_wise_cor_vis_aggxe)), 
                          rep("spot_wise", length(spot_wise_cor_vis_aggxe))))
ggplot(df, aes(x = type, y = cor)) + 
  geom_violin(trim=FALSE)


