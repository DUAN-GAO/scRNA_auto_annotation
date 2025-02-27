library(reticulate)

library(SingleCellExperiment)

library(SingleR)

library(scater)

ad = import("anndata")
sklearn <- import("sklearn")

train_file = "C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/Baron_human_test.h5ad"
test_file = "C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/Baron_human_train.h5ad"


train_ad = ad$read_h5ad(train_file)

counts_matrix = as.matrix(train_ad$X)

ref_sce = SingleCellExperiment::SingleCellExperiment(assays=list(counts=t(counts_matrix)))

rownames(ref_sce) = rownames(train_ad$var)


ref_sce = scater::logNormCounts(ref_sce)

ref_sce$Celltype = train_ad$obs$Celltype


test_ad = ad$read_h5ad(test_file)

test_counts_matrix = as.matrix(test_ad$X)

test_sce = SingleCellExperiment::SingleCellExperiment(assays=list(counts=t(test_counts_matrix)))

rownames(test_sce) = rownames(test_ad$var)

test_sce = scater::logNormCounts(test_sce)

pred <- SingleR(test = test_sce, ref = ref_sce, labels = ref_sce$Celltype)

acc = sklearn$metrics$accuracy_score(test_ad$obs$Celltype,pred$labels)

f1 = sklearn$metrics$f1_score(test_ad$obs$Celltype, pred$labels, average='macro')
