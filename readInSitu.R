load_as_spatial_experiment <- function(
    counts_file,
    metadata_file,
    metadata = NULL,
    id_columns = c("fov", "cell_ID"),
    coord_names = c("CenterX_global_px",
                    "CenterY_global_px"),
    ...
) {
  counts <- data.table::fread(counts_file)
  if (is.null(metadata))
    metadata <- data.table::fread(metadata_file)
  
  if (!is.null(id_columns)) {
    counts <- merge(metadata[, ..id_columns], counts)
    metadata <- merge(counts[, ..id_columns], metadata)
    counts <- counts[, -..id_columns]
  }
  
  spe <- SpatialExperiment::SpatialExperiment(
    assay = list(counts = t(counts)),
    colData = metadata,
    spatialCoordsNames = coord_names#,
    #...
  )
  
  colnames(counts(spe)) <- rownames(colData(spe)) <- rownames(spatialCoords(spe)) <- 1:ncol(counts(spe))
  
  return(spe)
}

load_cosmx_as_spatial_experiment <- function(
    counts_file,
    metadata_file,
    ...
) {
  load_as_spatial_experiment(
    counts_file,
    metadata_file,
    ...
  )
}

load_merscope_as_spatial_experiment <- function(
    counts_file,
    metadata_file,
    ...
) {
  load_as_spatial_experiment(
    counts_file,
    metadata = data.table::fread(metadata_file) |>
      setnames("V1", "cell"),
    id_columns = c("cell"),
    coord_names = c("center_x", "center_y"),
    ...
  )
}

load_xenium_as_spatial_experiment <- function(
    counts_file,
    metadata_file,
    ...
) {
  load_as_spatial_experiment(
    counts_file,
    metadata_file,
    id_columns = NULL,
    coord_names = c("x_centroid", "y_centroid"),
    ...
  )
}

.extract_indices <- function(idx, new.start, zero.based = TRUE) {
  if (length(idx) < 1) {
    return(NULL)
  }
  
  idx.cnts <- do.call(
    rbind,
    lapply(
      seq_len(length(new.start))[-1],
      function(x) c(x - ifelse(zero.based, 2, 1), new.start[[x]] - new.start[[x - 1]])
    )
  )
  colnames(idx.cnts) <- c("id", "n")
  
  return(
    list(
      i = idx,
      j = as.integer(tidyr::uncount(tibble::as_tibble(idx.cnts), n)[[1]]),
      new.start = new.start
    )
  )
}
