
# 安装R包
# install.packages("Seurat")

#加载
library(Seurat)
library(dplyr)

#包的版本
# sessionInfo()
# [1] dplyr_1.1.4        Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-0          

#工作目录
setwd("F:/4scRNA/1code/2Seurat")

#读入
pbmc = readRDS("1pbmc.rds")

#要把Seurat对象，当做一个数据库，包含多种数据
#其中，最重要的就是包含了多个表达矩阵和细胞注释信息（类似于临床信息）
#当需要用到表达矩阵或细胞注释信息时，会使用默认信息
#不想要使用默认信息，就必须指定或修改默认信息

#用@,$符号依次取下层，也可以[[]]

#assays 储存着表达矩阵
#counts 存储原始数据，是稀疏矩阵
#data 存储Normalize() 规范化的data
#scale.data 存储 ScaleData()缩放后的data
#SCT 储存SCT标准化之后的data

#meta.data中储存着细胞注释信息（类似于临床信息）

#active.assay 储存着默认的矩阵名

#active.ident 储存着默认的细胞注释信息（类似于临床信息）

dim(pbmc)
#获取表达矩阵
a = as.matrix(GetAssayData(object = pbmc@assays$RNA, layer = "counts")[1:20,1:20])
#或者用slot表示layer
b = as.matrix(GetAssayData(object = pbmc@assays$RNA, layer = "data")[1:20,1:20])
c = as.matrix(GetAssayData(object = pbmc@assays$SCT, layer = "counts")[1:20,1:20])
d = as.matrix(GetAssayData(object = pbmc@assays$SCT, layer = "data")[1:20,1:20])

#修改默认的表达矩阵
DefaultAssay(pbmc)
#DefaultAssay(pbmc) = "RNA"
e = as.matrix(GetAssayData(object = pbmc)[1:20,1:20])

table(d == e)


#获取细胞注释信息（类似于临床信息）
f = pbmc@meta.data

table(pbmc@meta.data$orig.ident)
colnames(pbmc@meta.data)

#修改默认的细胞注释信息
Idents(pbmc)
#Idents(pbmc) = "nCount_RNA"


#总之，要把Seurat对象，当做一个数据库，包含多种数据
#其中，最重要的就是包含了多个表达矩阵和细胞注释信息（类似于临床信息）
#当需要用到表达矩阵或细胞注释信息时，会使用默认信息
#不想要使用默认信息，就必须指定或修改默认信息

