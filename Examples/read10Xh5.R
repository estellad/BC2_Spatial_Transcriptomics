read10Xh5 <- function(dirname, fname = "filtered_feature_bc_matrix.h5") {
  spatial_dir <- file.path(dirname, "spatial")
  h5_file <- file.path(dirname, fname)
  
  if (!dir.exists(spatial_dir)) {
    stop("Spatial directory does not exist:\n  ", spatial_dir)
  }
  
  if (!file.exists(h5_file)) {
    stop("H5 file does not exist:\n  ", h5_file)
  }
  
  colData <- .read_spot_pos(spatial_dir)
  
  non.zero.indices <- .extract_indices(
    h5read(h5_file, "matrix/indices"),
    h5read(h5_file, "matrix/indptr")
  )
  
  rowData <- h5read(h5_file, "matrix/features/id")
  
  .counts <- sparseMatrix(
    i = non.zero.indices$i,
    j = non.zero.indices$j,
    x = h5read(h5_file, "matrix/data"),
    dims = h5read(h5_file, "matrix/shape"),
    dimnames = list(
      rowData,
      h5read(h5_file, "matrix/barcodes")
    ),
    index1 = FALSE
  )
  .counts <- .counts[, rownames(colData)]
  
  sce <- SingleCellExperiment(
    assays = list(
      counts = .counts
    ),
    rowData = rowData
    # ,
    # colData = colData
  )
  
  # Remove spots with no reads for all genes.
  sce <- sce[, Matrix::colSums(counts(sce)) > 0]
  
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  
  sce
}