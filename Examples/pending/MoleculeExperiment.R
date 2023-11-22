# TODO: does not work, write my own function for direct read in of SPE in SPE package

# BiocManager::install("MoleculeExperiment")
# Read in Xenium
library(MoleculeExperiment)
library(ggplot2)

path <- "~/Desktop/raw_data_must/Xenium_rawallexceptometiff/outs/"
me <- MoleculeExperiment::readXenium(path, keepCols = "essential")
me

# transform ME to SPE object
spe <- MoleculeExperiment::countMolecules(me)
# Error in getClass(Class, where = topenv(parent.frame())) : 
#   “gTMatrix” is not a defined class
spe

repoDirtest <- system.file("extdata", package = "MoleculeExperiment")
repoDirtest <- paste0(repoDirtest, "/xenium_V1_FF_Mouse_Brain")

metest <- readXenium(repoDirtest, keepCols = "essential")
metest

spetest <- countMolecules(metest)
spetest

## Xenium V1 FF mouse brain subset 
path <- "~/Desktop/raw_data_must/Xenium_V1_FF_MouseBrain_Subset/"
memousesmall <- MoleculeExperiment::readXenium(path, keepCols = "essential")
memousesmall

# transform ME to SPE object
spemousesmall <- MoleculeExperiment::countMolecules(memousesmall)
# Error in getClass(Class, where = topenv(parent.frame())) : 
#   “gTMatrix” is not a defined class


spemousesmall

















