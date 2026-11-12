# 安装 beachmat（确保编译成功）
# BiocManager::install("beachmat")
# 再安装 SingleR
# devtools::install_github("dviraran/SingleR@v2.4.0")
# devtools::install_github("dviraran/singleR")
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install("BiocNeighbors")
library(SingleR)
library(celldex) # 1.18.0 BiocManager::install("celldex")
library(Seurat) # 5.3.0
library(pheatmap) # 1.0.13

##下载注释数据库
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se

# 如需本地缓存数据库，可先保存再用 load() 读入；否则直接调用 celldex 函数下载
# save(hpca.se, file = "HumanPrimaryCellAtlas_hpca.se_human.RData")
bpe.se <- BlueprintEncodeData()

load("pbmc_tutorial.RData")
pbmc

meta=pbmc@meta.data #pbmc的meta文件，包含了seurat的聚类结果
head(meta)

plot1 <- DimPlot(pbmc, reduction = "umap", label = TRUE)

plot2<-DimPlot(pbmc, reduction = "tsne",
               label = TRUE)
plot1 + plot2

#进行singleR注释(两种方式)
pbmc_for_SingleR <- GetAssayData(pbmc, layer = "data") ##获取标准化矩阵
#pbmc_for_SingleR = as.SingleCellExperiment(pbmc)
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, 
                     ref = hpca.se, 
                     labels = hpca.se$label.main) #
pbmc.hesc

#seurat 和 singleR的table表
table(pbmc.hesc$labels,meta$seurat_clusters)
pbmc@meta.data$labels <- pbmc.hesc$labels

print(DimPlot(pbmc, group.by = c("seurat_clusters", "labels"),reduction = "umap"))

# 使用BP和HPCA两个数据库综合注释，使用list函数读入多个数据库
pbmc3 <- pbmc
pbmc3.hesc <- SingleR(test = pbmc_for_SingleR, ref = list(BP=bpe.se, HPCA=hpca.se), 
                      labels = list(bpe.se$label.main, hpca.se$label.main)) 
table(pbmc3.hesc$labels,meta$seurat_clusters)

pbmc3@meta.data$labels <-pbmc3.hesc$labels

print(DimPlot(pbmc3, group.by = c("seurat_clusters", "labels"),reduction = "umap"))

# 基于scores within cells
print(plotScoreHeatmap(pbmc.hesc))

# 基于 per-cell “deltas”诊断
plotDeltaDistribution(pbmc.hesc, ncol = 3)

# 与cluster结果比较
tab <- table(label = pbmc.hesc$labels,
             cluster = meta$seurat_clusters)

pheatmap(log10(tab + 10))

