# devtools::install_github("scverse/anndataR")
library(anndataR)
library(SingleCellExperiment)
library(BayesSpace)

vgg16_clustered_h5ad <- "~/Desktop/SampleData/AnnDataImage/VisPreprint/vgg16_clustered_kmeanlouvain.h5ad"
adata <- read_h5ad(vgg16_clustered_h5ad, to = "InMemoryAnnData")

sce <- adata$to_SingleCellExperiment()
sce

ggplot(data.frame(colData(sce)), 
       aes(x = sce$array_col, y = sce$array_row, color = factor(sce$louvain))) +
  geom_point() +
  labs(color = "cluster") +
  theme_bw() + scale_y_reverse() + scale_x_reverse()
