chrom_path <- paste0(root_path, "Chromium/")
Counts <- Seurat::Read10X_h5(paste0(chrom_path, "filtered_feature_bc_matrix.h5"))

# QC directly on annotated Chromium, which should be QCed already

chrom_full <- CreateSeuratObject(counts = Counts)
dim(chrom_full)

chrom_full <- preprocess_SEU(chrom_full)
table(chrom_full$libsize_drop) # 0
table(chrom_full$mito_drop) # 2027
chromfull_lowgenecount_drop <- unlist(chrom_full[["RNA"]][["lowgenecount_drop"]])
table(chromfull_lowgenecount_drop) # 4113 

chrom_full_qcd <- chrom_full[!chromfull_lowgenecount_drop, !chrom_full$mito_drop]
dim(chrom_full_qcd) # 13969 28338

dim(chrom_full)# 18082 30365

chrom_full_overlap <- chrom_full[, colnames(chrom_full) %in% overlap_cell_barcode]
chrom_full_overlap$percent.mt <- PercentageFeatureSet(chrom_full_overlap, pattern = "^MT-")

range(chrom_full_overlap$percent.mt)$ 0.00000 14.91619

hist(chrom_full_overlap$percent.mt)
MT_flag_overlap <- isOutlier(chrom_full_overlap$percent.mt, type = "higher")
table(MT_flag_overlap)

hist(chrom_full$percent.mt)
table(chrom_full$mito_drop)
range(chrom_full$percent.mt) # 0 - 79 really need to QC!
range(chrom_full_qcd$percent.mt) # 0 6.15
range(chrom$percent.mt) # 0 6.48
range(chrom_full_overlap$percent.mt) # 0 14

# Just relying on overlapping (10X QC), it gives still some cells with high mito percent. 0 - 14
# Directly QC on raw 30365 cells, returns range of mito percent to 0-6
# QC on overlapped with raw, which is what will be using in BC2 tutorial is 0-6. 

# Apparently 10X did not QC spots with high mito by OSCA book is.Outlier() function, but rather by some manual decided threshold. 
# Here we standardize the procedure by mapping in annotation, and perform QC on top of the overlapped object. 

# Plot cell-type percentage 
# Minimal theme + blue fill color
chrom_raw <- readRDS("./intermediate_data/chrom_raw.rds")
chrom_qcd <- readRDS("./intermediate_data/chrom_qcd.rds")

plot_ordered_bar <- function(seu, decreasing = TRUE, ylabel = "Percentage", 
                             title = "Cell-type Frequency in Annotated Chromium"){
  CD <- seu@meta.data
  cnt <- plyr::count(CD$Annotation)
  CD$Annotation <- factor(CD$Annotation, 
                          levels = cnt$x[order(cnt$freq, decreasing = decreasing)])
  
  p <- ggplot(data = CD, aes(x = Annotation)) + 
    # geom_bar(aes(y = (..count..)/sum(..count..))) +
    geom_bar() +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    ylab(ylabel) + ggtitle(title)
  
  p
}

p1 <- plot_ordered_bar(chrom_raw, decreasing = TRUE)
p2 <- plot_ordered_bar(chrom_qcd, decreasing = TRUE)

p1
p2

hist(chrom_full$percent.mt, breaks = 500, xlim = c(0,20))

hist(chrom_full_overlap$percent.mt, breaks = 500, xlim = c(0,20))

hist(chrom_full$nCount_RNA)

# New decision: Chromium QC with 10X default, but plot of discard flags, 
# because we are throwing away a lot of invasive tumor cell type if we use isOutlier on raw 30365 cells, or on overlapped object 27460, both case cut off at 6% 
# If we use the 10X QC threshold, then we keep the threshold up until 14, so it's much better
 

