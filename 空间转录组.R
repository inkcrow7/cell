rm(list = ls())
library(data.table)
library(RColorBrewer)
library(anndata) # reticulate 1.43.0
library(Matrix)
library(ggplot2)
library(future)
library(dplyr)
library(Seurat)
library(SeuratData) # devtools::install_github('satijalab/seurat-data')
library(cowplot)
library(patchwork)
library(jsonlite)
library(stringr)

# 读取MEX格式的数据
# filtered_feature_bc_matrix
# ├── barcodes.tsv.gz
# ├── features.tsv.gz
# └── matrix.mtx.gz
# 这3个文件必须是有固定的格式以及固定好的文件名字，不能修改！！！
library(data.table)
# 读取 matrix.mtx.gz 文件
matrix <- readMM("matrix.mtx.gz")
# 读取 barcodes.tsv.gz 文件
barcodes_df <- fread("barcodes.tsv.gz", header = FALSE, col.names = c("barcodes"))
# 读取 features.tsv.gz 文件
features_df <- fread("features.tsv.gz", header = FALSE, col.names = c("features"))
# 打印示例数据
cat("Matrix Dimensions:", dim(matrix), "\n")
cat("Barcodes DataFrame:\n")
head(barcodes_df)
cat("Features DataFrame:\n")
head(features_df)

## 先读取表达矩阵
# https://rdrr.io/github/satijalab/seurat/src/R/preprocessing.R
ct <- data.table::fread('GSE242249_RAW/GSM7757059_A4_matrix.csv.gz', data.table = F)
ct[1:4,1:4]
colnames(ct) <- gsub('[.]4','-1',colnames(ct))
rownames(ct) <- ct[,1]
ct <- ct[,-1]
ct[1:4,1:4]
dim(ct)

# 创建Seurat对象
spatial <- CreateSeuratObject(counts = ct, assay = "Spatial")
spatial
head(colnames(spatial))

## 读取HE切片信息
list.files('GSE242249_RAW/spatial/')
# https://satijalab.org/seurat/reference/read10x_image
# should include files tissue_lowres_iamge.png, scalefactors_json.json and tissue_positions_list.csv
# 注意tissue_positions_list.csv数据为6列，无表头
img <- Read10X_Image(image.dir = 'D:\\SG\\单细胞\\空间转录组\\GSE242249_RAW\\spatial',filter.matrix = TRUE)
img
head(img$centroids@coords)
tail(img$centroids@coords)

# 取交集
img <- img[Cells(x = spatial)]
DefaultAssay(spatial = img) <- "Spatial"
spatial[['slice1']] <- img
spatial

plot1 <- VlnPlot(spatial, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot1
plot2 <- SpatialFeaturePlot(spatial, features = "nCount_Spatial",pt.size.factor = 2.6) + theme(legend.position = "right")
plot2
pdf("wrap_plots.pdf",width = 7,height = 6)
wrap_plots(plot1, plot2)
dev.off()

## 标准化降为聚类分群（标准的降维聚类分群：空转数据推荐使用SCT标准化）
spatial <- SCTransform(spatial, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE )
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE, pcs.compute=50)
spatial <- FindNeighbors(spatial, reduction = "pca", dims = 1:30)
spatial <- FindClusters(spatial, verbose = FALSE, resolution = c(seq(0.1,1.5,by=0.1),2))
spatial <- RunUMAP(spatial, reduction = "pca", dims = 1:30)

# 基因表达可视化
SpatialFeaturePlot(spatial, features = c("PLCH2", "ACOT7"))
# pt.size.factor-这将扩展点的大小。默认值为1.6
# alpha-最低和最大的透明度。默认值为 c（1，1）。
# 尝试设置为 c（0.1， 1），以降低低表达点的透明度alpha
p1 <- SpatialFeaturePlot(spatial, features = "ACOT7", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(spatial, features = "ACOT7", alpha = c(0.1, 1))
p1 + p2

head(spatial@meta.data)
Idents(spatial) <- "SCT_snn_res.0.5"

p1 <- DimPlot(spatial, reduction = "umap", label = TRUE,pt.size = 2)
p1
p2 <- SpatialDimPlot(spatial,label = TRUE, label.size = 3, pt.size.factor = 2.3)
p2
pdf("SpatialDimPlot.pdf",width = 8,height = 7)
p1 + p2
dev.off()

SpatialDimPlot(spatial, cells.highlight = CellsByIdentities(object = spatial, 
                                                            idents = c(2, 1, 4, 3, 5, 8)), 
               facet.highlight = TRUE, ncol = 3)

# 空间高变基因的识别
de_markers <- FindMarkers(spatial, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = spatial, 
                   features = rownames(de_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)
