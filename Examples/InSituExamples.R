# library(BayesSpace)
library(R.utils)
library(rhdf5)
library(Matrix)
library(tibble)
library(tidyr)
library(Seurat) # CRAN package 
raw_data_path <- "~/Desktop/SampleData/Raw/"
tenxxenpath <- "10X_Xenium/Xenium/"

# Xe - Preprint original download
xeniumpath <- "~/Desktop/SampleData/archive/Xenium_Preprint_Data/Xenium/outs/"
# Xe - Mouse brain data original download 
xemousepath <- "/home/estelladong/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs/"

cosmx_lp9s1 <- "Nanostring_CosMx/Lung/P9/S1/"
merscope_path2 <- "Vizgen_MERSCOPE/Ovarian/P2/S2/"


## Xenium ------------------------------------------------------------
# count matrix
dirname <- xeniumpath
h5fname <- "cell_feature_matrix.h5"
coordfname <- "cells.csv.gz"

readXeniumSPE <- function(dirname, h5fname, coordfname, 
                          coord_names = c("x_centroid", "y_centroid")){
  h5_file <- file.path(dirname, h5fname)
  coord_file <- file.path(dirname, coordfname)
  
  # Count matrix 
  counts <- Seurat::Read10X_h5(h5_file)$`Gene Expression`
  
  # rowData
  featuresh5 <- h5read(h5_file, "matrix/features")
  rowData <- data.frame(gene_id = featuresh5$id[featuresh5$feature_type == "Gene Expression"],
                        gene_name = featuresh5$name[featuresh5$feature_type == "Gene Expression"])
  
  # Spatial and colData
  colData <- read.csv(gzfile(coord_file))
  
  spe <- SpatialExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData,
    spatialCoordsNames = coord_names
  )
  
  return(spe)
}


# cell coordinate names are "x_centroid" and "y_centroid" for Xenium
# feel free to rename after

# @Example Preprint
xe_spe <- readXeniumSPE(dirname = xeniumpath, 
                        h5fname = "cell_feature_matrix.h5", 
                        coordfname = "cells.csv.gz", 
                        coord_names = c("x_centroid", "y_centroid"))

# @Example MouseBrain
xe_spe <- readXeniumSPE(dirname = xemousepath, 
                        h5fname = "cell_feature_matrix.h5", 
                        coordfname = "cells.csv.gz", 
                        coord_names = c("x_centroid", "y_centroid"))




## CosMx -----------------------------------------------------------------
# count matrix
dirname <- "/home/estelladong/Téléchargements/Lung9_Rep1/Lung9_Rep1/Lung9_Rep1"
countmatfname <- "Lung9_Rep1_exprMat_file.csv"
metadatafname <- "Lung9_Rep1_metadata_file.csv"

readCosMxSPE <- function(dirname = dirname, 
                         countmatfname = countmatfname, 
                         metadatafname = metadatafname, 
                         coord_names = c("CenterX_global_px",
                                         "CenterY_global_px")){
  
  countmat_file <- file.path(dirname, countmatfname)
  metadata_file <- file.path(dirname, metadatafname)
  
  # Count matrix 
  countmat <- read.csv(countmat_file)
  countmat <- # TODO: harmonize names 
  
  # rowData
  
  # colData
  metadata <- read.csv(metadata_file)
  metadata <- 
  
  spe <- SpatialExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData,
    spatialCoordsNames = coord_names
  )
  
  return(spe)
}




to_align <- load_cosmx_as_spatial_experiment(
  counts_file = paste0(raw_data_path, cosmx_lp9s1, "cell_by_gene.csv"),
  metadata_file = paste0(raw_data_path, cosmx_lp9s1, "cell_metadata.csv")
)


# Merscope 
to_align2 <- load_merscope_as_spatial_experiment(
  counts_file = paste0(raw_data_path, merscope_path2, "OvarianP2S2_cell_by_gene.csv"),
  metadata_file = paste0(raw_data_path, merscope_path2, "OvarianP2S2_cell_metadata.csv")
)






library(BayesSpace)
sce <- exampleSCE()
clusterPlot(sce)



