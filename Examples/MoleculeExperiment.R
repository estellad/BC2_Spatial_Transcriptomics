library(MoleculeExperiment)
library(ggplot2)

path <- "~/Desktop/raw_data_must/Xenium/outs/"
me <- readXenium(path,
                 keepCols = "essential")
me
