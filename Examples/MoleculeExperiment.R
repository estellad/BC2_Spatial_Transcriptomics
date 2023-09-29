# BiocManager::install("MoleculeExperiment")
# Read in Xenium
library(MoleculeExperiment)
library(ggplot2)

path <- "~/Desktop/raw_data_must/Xenium_forMoleculeExperiment/outs/"
me <- MoleculeExperiment::readXenium(path, keepCols = "essential")
me

# transform ME to SPE object
spe <- MoleculeExperiment::countMolecules(me)
spe

repoDir <- system.file("extdata", package = "MoleculeExperiment")
repoDir <- paste0(repoDir, "/xenium_V1_FF_Mouse_Brain")

me <- readXenium(repoDir, keepCols = "essential")
me

spe <- countMolecules(me)
spe

















