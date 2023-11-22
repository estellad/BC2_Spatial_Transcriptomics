## SCE to SPE ------------------------------------------------------
spe <- as(sce, "SpatialExperiment")
spe@int_colData@listData$spatialCoords <- as.matrix(colData(sce)[, coord_names])

## SPE to SCE ------------------------------------------------------
colData(spe) <- cbind(colData(spe), spatialCoords(spe))
sce_new <- as(spe, "SingleCellExperiment")
