
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(dplyr)

#包的版本
# sessionInfo()
# [1] dplyr_1.1.4        patchwork_1.3.2    ggplot2_4.0.1      tibble_3.3.1      
# [5] stringr_1.6.0      data.table_1.18.0  Seurat_5.4.0       SeuratObject_5.3.0
# [9] sp_2.2-0      

#目录
setwd("F:/4scRNA/1code/6CellCycle")

#读取
pbmc = readRDS("F:/4scRNA/1code/5double/2pbmc_double.rds")

#细胞周期评分
pbmc <- NormalizeData(pbmc)

#获取G2M期相关基因
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(pbmc))

#获取S期相关基因
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(pbmc))

#细胞周期阶段评分
pbmc <- CellCycleScoring(pbmc, g2m.features=g2m_genes, s.features=s_genes)

colnames(pbmc@meta.data)
table(pbmc$Phase)

#画图
DimPlot(pbmc,group.by = "Phase",reduction = "tsne")

#保存
saveRDS(pbmc,"3pbmc_CellCycle.rds")

