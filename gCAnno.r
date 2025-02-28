library(reticulate)
library(openxlsx)

use_condaenv("gcANO")

build_model_script = "C:/Users/Administrator/Desktop/work/gCAnno/buildModelPipline.py"
define_celltype_bayes_script = "C:/Users/Administrator/Desktop/work/gCAnno/defineCellTypePipline_bayes.py"

# ad = import("anndata")
sklearn <- import("sklearn")

train_file = "C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/Zheng_68K_test.h5ad"
test_file = "C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/Zheng_68K_train.h5ad"

train_ad = ad$read_h5ad(train_file)
counts_matrix = t(as.matrix(train_ad$X))
colnames(counts_matrix) = gsub("-", "_", rownames(train_ad$obs))
rownames(counts_matrix) = rownames(train_ad$var)

#write.xlsx(counts_matrix, "expdata_normal.xls", rowNames = TRUE)
#write.csv(as.data.frame(counts_matrix), "expdata_normal.csv", rowNames = TRUE)


cell.name = rownames(train_ad$obs)
cluster = as.vector(train_ad$obs$Celltype)
dat = data.frame(cell.name,cluster)

#write.xls(dat, "cluster.xls")

output_dir <- "./result"
expdata_file <- "C:/Users/Administrator/Documents/expdata_normal.xls"
cluster_file <- "C:/Users/Administrator/Documents/cluster.xls"
feature_ratio <- 0.3

expdata_file <- "C:/Users/Administrator/Desktop/work/gCAnno/testData/test/expdata_normal.xls"
cluster_file <- "C:/Users/Administrator/Desktop/work/gCAnno/testData/test/cluster.xls"

command_build_model <- paste(
    "python", build_model_script,
    "-o", output_dir,
    "-e", expdata_file,
    "-c", cluster_file,
    "-fr", feature_ratio
)

system(command_build_model)
