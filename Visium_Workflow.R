# Visium
library(scater)
sce <- mockSCE()

is_mito <- grepl("(^MT-)|(^mt-)", rowData(sce)$symbol)
table(is_mito)

# calculate per-spot QC metrics and store in colData
sce <- scuttle::addPerCellQC(sce)
