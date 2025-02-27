library(SingleCellExperiment)
library(CHETAH)
library(reticulate)
ad = import("anndata")
sklearn <- import("sklearn")

train_file = "C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/Zheng_68K_test.h5ad"
test_file = "C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/Zheng_68K_train.h5ad"

train_ad = ad$read_h5ad(train_file)
counts_matrix = t(as.matrix(train_ad$X))
colnames(counts_matrix) = rownames(train_ad$obs)
rownames(counts_matrix) = rownames(train_ad$var)
ref_ct = as.vector(train_ad$obs$Celltype)
reference <- SingleCellExperiment(
  assays = list(counts = counts_matrix),
  colData = data.frame(celltypes = ref_ct)
)

test_ad = ad$read_h5ad(test_file)
test_counts_matrix = t(as.matrix(test_ad$X))
colnames(test_counts_matrix) = rownames(test_ad$obs)
rownames(test_counts_matrix) = rownames(test_ad$var)
ref_ct = as.vector(test_ad$obs$Celltype)
input <- SingleCellExperiment(
  assays = list(counts = test_counts_matrix)
)

result <- CHETAHclassifier(input = input,
                              ref_cells = reference)
pred = as.vector(result$celltype_CHETAH)

acc = sklearn$metrics$accuracy_score(test_ad$obs$Celltype,pred)
f1 = sklearn$metrics$f1_score(test_ad$obs$Celltype, pred, average='macro')