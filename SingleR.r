library(reticulate)

library(SingleCellExperiment)

library(SingleR)

library(scater)

train_file = "C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/Zheng_68K_test.h5ad"

train_ad = ad$read_h5ad(train_file)

counts_matrix = as.matrix(train_ad$X)

ref_sce=SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts_matrix))

colData(ref_sce)$Type=colnames(train_ad$obs$Celltype)

ref_sce = scater::logNormCounts(ref_sce)