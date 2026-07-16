
# 安装R包
# install.packages("devtools")
# devtools::install_github("PaulingLiu/ROGUE")
# install.packages("clustree")
# install.packages("harmony")

#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ROGUE)
library(clustree)
library(harmony)
library(dplyr)

#包的版本
# sessionInfo()
# [1] dplyr_1.1.4        harmony_1.2.4      Rcpp_1.1.1         clustree_0.5.1    
# [5] ggraph_2.2.2       ROGUE_1.0          patchwork_1.3.2    ggplot2_4.0.1     
# [9] tibble_3.3.1       stringr_1.6.0      data.table_1.18.0  Seurat_5.4.0      
# [13] SeuratObject_5.3.0 sp_2.2-0      

#目录
setwd("F:/4scRNA/1code/8harmony")

#读取
pbmc = readRDS("F:/4scRNA/1code/7Normalize/4pbmc_Normalize.rds")

#设置激活的矩阵及分组信息
DefaultAssay(pbmc)
DefaultAssay(pbmc) = "SCT"
DefaultAssay(pbmc)

table(Idents(pbmc))
Idents(pbmc) = "orig.ident"
table(Idents(pbmc))

#10个变化最大的基因
top10 <- head(VariableFeatures(pbmc), 10)

#画图
pdf(file="1.pdf",width=7,height=6)
VariableFeaturePlot(object = pbmc)
dev.off()

pdf(file="2.pdf",width=7,height=6)
LabelPoints(plot = VariableFeaturePlot(object = pbmc), points = top10, repel = TRUE)
dev.off()

#PCA
pbmc <- RunPCA(pbmc, verbose = F)

#主成分分析图形
pdf(file="3.pdf",width=7,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#绘制每个PCA成分的相关基因
pdf(file="4.pdf",width=10,height=9)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#主成分分析热图
pdf(file="5.pdf",width=10,height=9)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#选取合适的PC
#主成分累积贡献大于90%,选择拐点
pdf(file="6.pdf",width=7,height=6)
ElbowPlot(pbmc, ndims = 50)
dev.off()

#确定与每个 PC 的百分比   
pct <- pbmc [["pca"]]@stdev / sum( pbmc [["pca"]]@stdev) * 100
pct

#计算每个 PC 的累计百分比
cumu <- cumsum(pct)
cumu

#设置PC
pcs = 1:40

#harmony
pbmc <- RunHarmony(pbmc, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)

table(pbmc@meta.data$orig.ident)

#选取合适的分辨率
#从0.1-2的resolution结果均运行一遍
seq = seq(0.1,2,by=0.1)
pbmc <- FindNeighbors(pbmc,  dims = pcs) 
for (res in seq){
  pbmc = FindClusters(pbmc, resolution = res)
}

#画图
p1 = clustree(pbmc,prefix = "SCT_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
ggsave("SCT_sun_res.png", p, width = 30, height = 14)

#降维聚类
# pbmc <- FindNeighbors(pbmc, reduction = "harmony",  dims = pcs) %>% FindClusters(resolution = 1)
# pbmc <- RunUMAP(pbmc, reduction = "harmony",  dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
#降维聚类
pbmc <- FindNeighbors(pbmc, reduction = "pca",  dims = pcs) %>% FindClusters(resolution = 1)
pbmc <- RunUMAP(pbmc, reduction = "pca",  dims = pcs) %>% RunTSNE(dims = pcs, reduction = "pca")

colnames(pbmc@meta.data)
#画图
pdf(file="7.pdf",width=7,height=6)
DimPlot(pbmc, reduction = "umap", label = T)
dev.off()

pdf(file="8.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "umap",label = F,group.by = "orig.ident")
dev.off()

pdf(file="9.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "umap",label = F,group.by = "Is_Double")
dev.off()

pdf(file="10.pdf",width=7,height=6)
DimPlot(pbmc, reduction = "tsne", label = T)
dev.off()

pdf(file="11.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "tsne",label = F,group.by = "orig.ident")
dev.off()

pdf(file="12.pdf",width=7,height=6)
DimPlot(pbmc,reduction = "tsne",label = F,group.by = "Is_Double")
dev.off()

#保存
saveRDS(pbmc,"5pbmc_UMPA.TSNE.rds")

