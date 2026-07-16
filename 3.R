
# 安装R包
# install.packages("data.table")
# install.packages("stringr")
# install.packages("tibble")

#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(dplyr)

#包的版本
# sessionInfo()
# [1] tibble_3.3.1       stringr_1.6.0      data.table_1.18.0  dplyr_1.1.4       
# [5] Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-0    

#工作目录
setwd("F:/4scRNA/1code/3Read")


#####1.matrix.mtx、genes.tsv和barcodes.tsv####

#目录
setwd("F:/4scRNA/1code/3Read/1GSE234527")

# 获取数据文件夹下的所有样本文件列表
samples <- list.files("seurat/")
samples

# 创建一个空的列表
seurat_list <- list()

#读取数据并创建Seurat对象
#删除，小于200个基因表达的细胞，小于3个细胞表达的基因
for (sample in samples) {
  #文件路径
  data.path <- paste0("seurat/", sample)
  
  #读取10x数据
  seurat_data <- Read10X(data.dir = data.path)

  #创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = seurat_data,project = sample,min.features = 200,min.cells = 3)
  
  #添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}

#合并
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)

#将Layers融合
pbmc = JoinLayers(seurat_combined)

#####2.H5格式#####

#目录
setwd("F:/4scRNA/1code/3Read/2GSE199866")

#获取文件名称
fs=list.files(pattern = '.h5')
fs

#依次读入并创建Seurat对象
sceList = lapply(fs, function(x){
  a=Read10X_h5(x)
  p=str_split(x,'_',simplify = T)[,1]
  sce <- CreateSeuratObject(a,project = p ,min.features = 200,min.cells = 3)
})

#获取样本GSM号
folders = substr(fs,1,10)
folders

#使用merge函数进行合并
sce.big <- merge(sceList[[1]], 
                y = sceList[-1], 
                add.cell.ids = folders)

#将Layers融合
pbmc = JoinLayers(sce.big)

#####3.R数据文件(RDS/RDATA文件)####

#目录
setwd("F:/4scRNA/1code/3Read/3")

#读取RDATA文件
load(file ="1pbmc.RData")

#读取RDS文件
pbmc2 = readRDS("1pbmc.rds")


#####4.TXT或CSV####

#目录
setwd("F:/4scRNA/1code/3Read/4/GSE153935")

#读入，用fread，由于矩阵很大，而且格式很混乱
#一定要用fread，否则其他函数要很久
pbmc <-fread("GSE153935_TLDS_AllCells.txt", sep="\t")
#pbmc[1:5,1:5]

#直接转化为行名
pbmc <-  column_to_rownames(pbmc,"V1")
#pbmc[1:5,1:5]

#创建seurat对象
pbmc <- CreateSeuratObject(pbmc, min.features = 300, min.cells = 3)

#将Layers融合
pbmc = JoinLayers(pbmc)





#目录
setwd("F:/4scRNA/1code/3Read/4/GSE165722")

#获取样本位置
fs2=list.files(pattern = '.txt')
fs2
fs3=list.files(pattern = '.tsv')
fs3

#获取样本名称
folders = substr(fs2,1,10)
folders

#读入
sceList = list()
for (i in 1:length(fs2)) {
  #读入
  abc123 = as.data.frame(fread(fs2[i]))
  abc456 = as.data.frame(fread(fs3[i]))
  #将gene列转化为行名，基因名
  abc456 <-  column_to_rownames(abc456,"gene")
  #添加矩阵的列名，细胞名
  colnames(abc456) = abc123[,1]
  #创建seurat对象
  sce <- CreateSeuratObject(abc456,project = folders[i],min.features = 300, min.cells = 3)
  #依次放入list中
  sceList[i] = sce
}


#合并
sce.big <- merge(sceList[[1]], 
                 y = sceList[-1], 
                 add.cell.ids = folders)

#将Layers融合
pbmc = JoinLayers(sce.big)



#添加细胞注释信息

#获取细胞注释信息
abc789 = pbmc@meta.data

#导出
write.table(data.frame(ID=rownames(abc789),abc789),file="meta.txt", sep="\t", quote=F, row.names = F,col.names = T)

#读入
meta = fread("meta.xlsx")

#将第一列转化为行名
meta <-  column_to_rownames(meta,"ID")

#添加meta.data信息，细胞名的顺序必须一致
pbmc <- AddMetaData(object = pbmc, 
                    metadata = meta,   
                    col.name = c("group1","Type2")) 

#当不一致时，根据seurat对象的名字进行排序
meta = meta[colnames(pbmc),]

