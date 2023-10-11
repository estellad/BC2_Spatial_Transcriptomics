BiocManager::install("STdeconvolve")
library(STdeconvolve)

data(mOB)
pos <- mOB$pos ## x and y positions of each pixel
cd <- mOB$counts ## matrix of gene counts in each pixel
annot <- mOB$annot ## annotated tissue layers assigned to each pixel

# STdeconvolve first feature selects for genes most likely to be relevant for 
# distinguishing between cell-types by looking for highly overdispersed genes across ST pixels. 
# Pixels with too few genes or genes with too few reads can also be removed.

## remove pixels with too few genes
counts <- cleanCounts(counts = cd,
                      min.lib.size = 100,
                      min.reads = 1,
                      min.detected = 1,
                      verbose = TRUE)

## feature select for genes
corpus <- restrictCorpus(counts,
                         removeAbove = 1.0,
                         removeBelow = 0.05,
                         alpha = 0.05,
                         plot = TRUE,
                         verbose = TRUE)

## Note: the input corpus needs to be an integer count matrix of pixels x genes
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)

## select model with minimum perplexity
optLDA <- optimalModel(models = ldas, opt = "min")

results <- getBetaTheta(optLDA,
                        perc.filt = 0.05,
                        betaScale = 1000)

deconProp <- results$theta
deconGexp <- results$beta

## Visualization ----------------------------------------------------
vizAllTopics(deconProp, pos, 
             #groups = annot, 
             #group_cols = rainbow(length(levels(annot))),
             r=0.4)

vizAllTopics(deconProp, pos, 
             groups = annot, 
             group_cols = rainbow(length(levels(annot))),
             r=0.4)

# For faster plotting, we can visualize the pixel proportions of a single 
# cell-type separately using vizTopic()
vizTopic(theta = deconProp, pos = pos, topic = "8", plotTitle = "X8",
         size = 5, stroke = 1, alpha = 0.5,
         low = "white",
         high = "red")

## Still want to infer transcriptional profile -------------------------
# proxy theta for the annotated layers
mobProxyTheta <- model.matrix(~ 0 + annot)
rownames(mobProxyTheta) <- names(annot)
# fix names
colnames(mobProxyTheta) <- unlist(lapply(colnames(mobProxyTheta), function(x) {
  unlist(strsplit(x, "annot"))[2]
}))

mobProxyGexp <- counts %*% mobProxyTheta


# 2.3 Annotation Strategy 1: Transcriptional correlations
corMtx_beta <- getCorrMtx(# the deconvolved cell-type `beta` (celltypes x genes)
  m1 = as.matrix(deconGexp),
  # the reference `beta` (celltypes x genes)
  m2 = t(as.matrix(mobProxyGexp)),
  # "b" = comparing beta matrices, "t" for thetas
  type = "b")


## row and column names need to be characters
rownames(corMtx_beta) <- paste0("decon_", seq(nrow(corMtx_beta)))

correlationPlot(mat = corMtx_beta,
                # colLabs (aka x-axis, and rows of matrix)
                colLabs = "Deconvolved cell-types",
                # rowLabs (aka y-axis, and columns of matrix)
                rowLabs = "Ground truth cell-types",
                title = "Transcriptional correlation", annotation = TRUE) +
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))


# ---------------------------------------------------
corMtx_theta <- getCorrMtx(# deconvolved cell-type `theta` (pixels x celltypes)
  m1 = as.matrix(deconProp),
  # the reference `theta` (pixels x celltypes)
  m2 = as.matrix(mobProxyTheta),
  # "b" = comparing beta matrices, "t" for thetas
  type = "t")


## row and column names need to be characters
rownames(corMtx_theta) <- paste0("decon_", seq(nrow(corMtx_theta)))

correlationPlot(mat = corMtx_theta,
                # colLabs (aka x-axis, and rows of matrix)
                colLabs = "Deconvolved cell-types",
                # rowLabs (aka y-axis, and columns of matrix)
                rowLabs = "Ground truth cell-types",
                title = "Proportional correlation", annotation = TRUE) +
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))


## Order the cell-types rows based on best match (highest correlation) with
## each community.
## Cannot have more rows than columns for this pairing, so transpose
pairs <- lsatPairs(t(corMtx_theta))
m <- t(corMtx_theta)[pairs$rowix, pairs$colsix]

correlationPlot(mat = t(m), # transpose back
                # colLabs (aka x-axis, and rows of matrix)
                colLabs = "Deconvolved cell-types",
                # rowLabs (aka y-axis, and columns of matrix)
                rowLabs = "Ground truth cell-types",
                title = "Transcriptional correlation", annotation = TRUE) +
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))


# 2.4 Annotation Strategy 2: GSEA -----------------------------------------
mobProxyLayerMarkers <- list()

## make the tissue layers the rows and genes the columns
gexp <- t(as.matrix(mobProxyGexp))

for (i in seq(length(rownames(gexp)))){
  celltype <- i
  ## log2FC relative to other cell-types
  ## highly expressed in cell-type of interest
  highgexp <- names(which(gexp[celltype,] > 10))
  ## high log2(fold-change) compared to other deconvolved cell-types and limit
  ## to the top 200
  log2fc <- sort(
    log2(gexp[celltype,highgexp]/colMeans(gexp[-celltype,highgexp])),
    decreasing=TRUE)[1:200]
  
  ## for gene set of the ground truth cell-type, get the genes
  ## with log2FC > 1 (so FC > 2 over the mean exp of the other cell-types)
  markers <- names(log2fc[log2fc > 1])
  mobProxyLayerMarkers[[ rownames(gexp)[celltype] ]] <- markers
}

celltype_annotations <- annotateCellTypesGSEA(beta = results$beta,
                                              gset = mobProxyLayerMarkers,
                                              qval = 0.05)
celltype_annotations$results$`2`
celltype_annotations$predictions


####################################################################
#                             SPE Input                            #
####################################################################
# BiocManager::install("TENxVisiumData")
library(SpatialExperiment)
library(TENxVisiumData)

## load the MouseBrainCoronal SpatialExperiment object from `TENxVisiumData`
se <- TENxVisiumData::MouseBrainCoronal()

# (Alternatively, download the data from website and create a new directory
# First, make a directory to store the downloaded files:
  
f <- "visiumTutorial/"

if(!file.exists(f)){
  dir.create(f)
}

# Download and unzip the Feature / barcode matrix (filtered) and the Spatial imaging data:
if(!file.exists(paste0(f, "V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"))){
  tar_gz_file <- "http://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"
  download.file(tar_gz_file, 
                destfile = paste0(f, "V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"), 
                method = "auto")
}
untar(tarfile = paste0(f, "V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.tar.gz"), 
      exdir = f)

if(!file.exists(paste0(f, "V1_Adult_Mouse_Brain_spatial.tar.gz"))){
  spatial_imaging_data <- "http://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_spatial.tar.gz"
  download.file(spatial_imaging_data, 
                destfile = paste0(f, "V1_Adult_Mouse_Brain_spatial.tar.gz"), 
                method = "auto")
}
untar(tarfile = paste0(f, "V1_Adult_Mouse_Brain_spatial.tar.gz"), 
      exdir = f)
  

se <- SpatialExperiment::read10xVisium(samples = f,
                                       type = "sparse",
                                       data = "filtered")

















