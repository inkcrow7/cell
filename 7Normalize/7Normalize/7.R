
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(dplyr)
options(future.globals.maxSize = 1e12)

#包的版本
# sessionInfo()
# [1] dplyr_1.1.4        patchwork_1.3.2    ggplot2_4.0.1      tibble_3.3.1      
# [5] stringr_1.6.0      data.table_1.18.0  Seurat_5.4.0       SeuratObject_5.3.0
# [9] sp_2.2-0   

#目录
setwd("F:/4scRNA/1code/7Normalize")

#读取
pbmc = readRDS("F:/4scRNA/1code/6CellCycle/3pbmc_CellCycle.rds")

# 每个细胞最初包含相同数量的 RNA 分子的假设。
# 结果存储在pbmc@assays$RNA@layers$data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 由于单细胞是稀疏矩阵，很多基因的表达值几乎为0
# 有助于突出单细胞数据集中的生物信号
# 寻找2000个高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 改变每个基因的表达，使细胞间的平均表达为 0
# 缩放每个基因的表达，使细胞间的方差为 1
# 在下游分析中给予同等权重。
# 在这里去除细胞周期的影响
# 结果存储在pbmc@assays$RNA@layers$scale.data
pbmc <- ScaleData(pbmc,vars.to.regress = c("S.Score", "G2M.Score"))
#pbmc <- ScaleData(pbmc,vars.to.regress = c("S.Score", "G2M.Score"),features = rownames(pbmc))

# SCT，相当于替代了上述的三个函数NormalizeData，FindVariable，ScaleData
# 寻找3000个高变基因

# 在常规分析中，使用少量的PC既能关注到关键的生物学差异，
# 又能够不引入更多的技术差异，相当于一种保守性的做法。
# 它会失去一些生物差异信息但是同时又在常规手段中比较安全。
# 但SCT的归一化、标准化都做得不错，
# 多输入一些PCs能提取更多的生物差异，并且兼顾不引入技术误差。
# SCT认为:新增加的这1000个基因就包含了之前没有检测到的微弱的生物学差异。
# 而且，即使使用全部的全部的基因去做下游分析，得到的结果也是和SCT的结果类似

# 结果存储在pbmc@assays$SCT中
pbmc <- SCTransform(pbmc, vars.to.regress = c("S.Score", "G2M.Score"))

#默认的矩阵
DefaultAssay(pbmc)
DefaultAssay(pbmc) = "RNA"
DefaultAssay(pbmc)

#后续还是以SCT为例进行
DefaultAssay(pbmc) = "SCT"
DefaultAssay(pbmc)

#保存
saveRDS(pbmc,"4pbmc_Normalize.rds")
