library(reticulate)

# train_file = "C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/Baron_human_test.h5ad"
# test_file = "C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/Baron_human_train.h5ad"
# use_condaenv("gcANO")

py_run_string("import anndata")
py_run_string("import pandas")

file <- "Zheng_68K" #唯一需要改变的文件名称变量
train_file_path <- paste0("C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/", file, "_test.h5ad")
test_file_path <- paste0("C:/Users/Administrator/Desktop/scRNA_annotation/SCSA/data/", file, "_train.h5ad")

# 在 Python 中设置文件路径
py_run_string(paste0("train_file = '", train_file_path, "'"))
py_run_string(paste0("test_file = '", test_file_path, "'"))

# 或者直接在 R 中构建路径，然后传递给 Python
py_run_string(paste0("train_file = '", train_file_path, "'"))
py_run_string(paste0("test_file = '", test_file_path, "'"))

# 使用pandas构建ref文件，训练模型时使用
py_run_string("train_ad = anndata.read_h5ad(train_file)")
py_run_string("counts_matrix = pandas.DataFrame(train_ad.X.transpose())")
py_run_string("counts_matrix.columns = train_ad.obs.index")
py_run_string("counts_matrix.index = train_ad.var.index")
py_run_string("counts_matrix.to_csv('C:/Users/Administrator/Documents/expdata_normal.csv', sep='\t', encoding='utf-8')")

# 使用pandas构建test文件，测试模型时使用
py_run_string("test_ad = anndata.read_h5ad(test_file)")
py_run_string("counts_matrix = pandas.DataFrame(test_ad.X.transpose())")
py_run_string("counts_matrix.columns = test_ad.obs.index")
py_run_string("counts_matrix.index = test_ad.var.index")
py_run_string("counts_matrix.to_csv('C:/Users/Administrator/Documents/test_normal.csv', sep='\t', encoding='utf-8')")

# 使用pandas构建cluster文件，训练模型时使用
py_run_string("cellname = train_ad.obs.index")
py_run_string("cluster = [s.replace(' ', '.') for s in list(train_ad.obs.Celltype)]")
py_run_string("cluster_data = {
    'cell.name': cellname,
    'cluster': cluster
}")
py_run_string("df = pandas.DataFrame(cluster_data)")
py_run_string("df.to_csv('C:/Users/Administrator/Documents/cluster.csv', sep='\t', encoding='utf-8', index=False)")




# 切换conda环境运行，注意这里可能会报错，需要重新启动R运行
library(reticulate)

use_condaenv("gcANO")
build_model_script = "C:/Users/Administrator/Desktop/work/gCAnno/buildModelPipline.py"

output_dir <- "C:/Users/Administrator/Documents/"
expdata_file <- "C:/Users/Administrator/Documents/expdata_normal.csv"
cluster_file <- "C:/Users/Administrator/Documents/cluster.csv"
feature_ratio <- 0.3

command_build_model <- paste(
    "python", build_model_script,
    "-o", output_dir,
    "-e", expdata_file,
    "-c", cluster_file,
    "-fr", feature_ratio
)

system(command_build_model)


define_celltype_bayes_script <- "C:/Users/Administrator/Desktop/work/gCAnno/defineCellTypePipline_bayes.py"

output_dir_bayes <- "C:/Users/Administrator/Documents/"
embedding_distance_file <- "C:/Users/Administrator/Documents/embedding/02.embeddingDis.xls"
query_expdata_file <- "C:/Users/Administrator/Documents/test_normal.csv"
sg_value <- 65

# 构建命令
command_annotate_bayes <- paste(
    "python", define_celltype_bayes_script,
    "-o", output_dir_bayes,
    "-d", embedding_distance_file,
    "-q", query_expdata_file,
    "-sg", sg_value
)

# 执行命令
system(command_annotate_bayes)

# 计算accuracy和f1评分
sklearn <- import("sklearn")

df = read.csv("out_predict_cluster.xls",sep = "")
ad = import("anndata")
test_ad = ad$read_h5ad(test_file_path)

acc = sklearn$metrics$accuracy_score(test_ad$obs$Celltype,df$pridect)
f1 = sklearn$metrics$f1_score(test_ad$obs$Celltype, df$pridect, average='macro')