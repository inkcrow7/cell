#devtools::install_github("sqjin/CellChat")
#devtools::install_github('satijalab/seurat-data')
#install.packages("Matrix")
#install.packages(
  #"https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.1.tar.gz",
  #repos = NULL,
  #type = "source"
#)
#install.packages("SeuratObject")
#install.packages("Seurat") 
library(Seurat)
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork) # 最强大的拼图包
library(ggplot2)
library(svglite)
options(stringsAsFactors = FALSE)

pbmc3k.final <- readRDS("pbmc.rds")
pbmc3k.final@commands$FindClusters  # 也可以看一看作者的其他命令，Seurat是记录其分析过程的。

cellchat <- createCellChat(object=pbmc3k.final,group.by = "cell_type")
#cellchat <- createCellChat(pbmc3k.final@assays$RNA@data, meta = pbmc3k.final@meta.data, group.by = "cell_type")
cellchat
summary(cellchat)
str(cellchat)
levels(cellchat@idents)
#cellchat <- setIdent(cellchat, ident.use = "cell_type")
groupSize <- as.numeric(table(cellchat@idents))  
#查看每个cluster有多少个细胞，后面画图的时候需要用到这个值
groupSize
# [1] 684 481 476 344 291 162 155  32  13

CellChatDB <- CellChatDB.human
#导入小鼠是CellChatDB <- CellChatDB.mouse
str(CellChatDB) #查看数据库信息
#包含interaction、complex、cofactor和geneInfo这4个dataframe
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new() #下一步不出图的时候运行
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种；提取向量或数据框列中的唯一值（去重后的值）
#选择"Secreted Signaling"进行后续细胞互作分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use # set the used database in the object

# 预处理
## 在矩阵的所有的基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个）
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#上一步运行的结果储存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human) 
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project

# 推断细胞通讯网络
# 根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")

# 通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")

# 细胞互作关系展示
# 7.1 所有细胞群总体观：细胞互作数量与强度统计分析：
# 统计细胞和细胞之间通信的数量（有多少个配体-受体对）和强度（概率）
cellchat <- aggregateNet(cellchat)
#计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=FALSE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
# save as TIL/net_number_strength.pdf

# 检查每种细胞发出的信号
#mat <- cellchat@net$count
#par(mfrow = c(3,3), xpd=TRUE)
par("mar")
par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as TIL/net_number_individual.pdf

## 运行上述代码出现报错 Error in plot.new() : figure margins too large
# par("mar")
## [1] 5.1 4.1 4.1 2.1
# par(mar=c(1,1,1,1))
# 重新运行上面的代码

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as TIL/net_strength_individual.pdf

mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as TIL/net_strength_individual.pdf

# 7.2 单个信号通路或配体-受体介导的细胞互作可视化（层次图、网络图、和弦图、热图）
# 比如在前面的功能富集分析或case control的比较中找到了一些信号通路差异，就可以进一步在细胞水平上验证。
cellchat@netP$pathways  #查看都有哪些信号通路
# 选择其中一个信号通路，比如说TGFb
pathways.show <- c("TGFb")  

# 层次图（Hierarchy plot）
levels(cellchat@idents)    # show all celltype
# [1] "Naive CD4 T"  "Memory CD4 T" "CD14+ Mono"   "B"            "CD8 T"       
# [6] "FCGR3A+ Mono" "NK"           "DC"           "Platelet"    
vertex.receiver = c(1,2,4,6) # 定义一个数字向量（淋系细胞），将细胞类型的索引作为目标
#par(mar=c(5.1,4.1,4.1,2.1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "hierarchy", vertex.receiver = vertex.receiver)
# save as TIL/CXCL_hierarchy.pdf；一定要定义layout = "hierarchy"。

# 网络图（Circle plot）
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# save as TIL/CXCL_circle.pdf

#和弦图（Chord diagram）
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord") # 需将窗口拉大
# save as TIL/CXCL_chord.pdf

# 热图（Heatmap）
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# save as TIL/CXCL_heatmap.pdf

# 7.3 配体-受体层级的可视化（计算各个ligand-receptor pair对信号通路的贡献）
#计算配体受体对选定信号通路的贡献值（在这里就是查看哪条信号通路对TGFb贡献最大）
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE) #提取对TGFb有贡献的所有配体受体 
# save as TIL/CXCL_LR_contribution.pdf

# 层次图
#提取对这个通路贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
LR.show <- pairLR.TGFb[1,] 
vertex.receiver = c(1,2,4,6) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show, layout = "hierarchy", pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# save as TIL/CXCL_hierarchy2.pdf

# 网络图（Circle plot）
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# save as TIL/CXCL_circle2.pdf

# 和弦图（Chord diagram）
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord") # 窗口拉大
# save as TIL/CXCL_chord2.pdf

# 7.4 自动（批量）保存每个信号通路的互作结果
# Access all the signaling pathways showing significant communications将所有信号通路找出来
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = c(1,2,4,6) #不画层次图就不需要这一步
dir.create("all_pathways_com_circle") #创建文件夹保存批量画图结果
setwd("all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"),
            vertex.receiver = vertex.receiver, layout = "circle") #绘制网络图
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
         plot=gg, width = 5, height = 2.5, units = 'in', dpi = 300)
}
setwd("../")

# 7.5多个配体-受体介导的细胞互作关系可视化
# 气泡图（全部配体受体）
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
#需要指定受体细胞和配体细胞
p = netVisual_bubble(cellchat, sources.use = c(3,5,7,8,9), 
                     targets.use = c(1,2,4,6), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble.pdf", p, width = 8, height = 12) #髓系对淋巴的调节
# save as TIL/Mye_Lymph_bubble.pdf

# 气泡图（指定信号通路或配体-受体）
#比如制定CCL和CXCL这两个信号通路
netVisual_bubble(cellchat, sources.use = c(3,5,7,8,9), targets.use = c(1,2,4,6), 
                 signaling = c("CCL","CXCL"), remove.isolate = FALSE)

# 气泡图（指定信号通路或配体-受体并指定细胞）
# show all the significant interactions (L-R pairs) based on user's input
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","TGFb"))
netVisual_bubble(cellchat, sources.use = c(3,6,8), targets.use = c(1,4,5), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)

# 参与某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示（小提琴图和气泡图）
## Plot the signaling gene expression distribution
p = plotGeneExpression(cellchat, signaling = "TGFb")
ggsave("TGFb_GeneExpression_vln.pdf", p, width = 8, height = 8)
p = plotGeneExpression(cellchat, signaling = "TGFb", type = "dot")
ggsave("TGFb_GeneExpression_dot.pdf", p, width = 8, height = 6)

# 通讯网络系统分析（可信度有待考量）
# 1 社会网络分析（通讯网络中的角色识别）
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                   width = 15, height = 6, font.size = 10)
# 识别细胞的信号流模式
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 5)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 5)
ht1 + ht2

# 非负矩阵分解（NMF）识别细胞的通讯模式（这里是一个比较标准的NMF的应用方式）
# 信号输出细胞的模式识别
#计算分解成几个因子(pattern)比较合适（这一步运行比较慢 。在使用NMF对细胞进行亚群细分时，如果不测试的话，最好选择比细胞类型多一点的值）
selectK(cellchat, pattern = "outgoing")

# 热图
nPatterns = 2 # 挑选曲线中第一个出现下降的点（从2就开始下降了）
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, 
                                          width = 5, height = 9, font.size = 6)

# 冲击图/河流图（river plot）
netAnalysis_river(cellchat, pattern = "outgoing")

# 气泡图
netAnalysis_dot(cellchat, pattern = "outgoing")

# 信号输入细胞的模式识别
selectK(cellchat, pattern = "incoming") 

# 热图
nPatterns = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, 
                                          width = 5, height = 9, font.size = 6)

# 冲击图
netAnalysis_river(cellchat, pattern = "incoming")

# 气泡图
netAnalysis_dot(cellchat, pattern = "incoming")

# 信号网络的流行学习与分类
# 把共同起作用的信号通路归纳在一起，分为基于功能的归纳和基于拓扑结构的归纳
# 基于功能相似性的流行学习与分类
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#Error in do_one(nmeth) : NA/NaN/Inf in foreign function call (arg 1)
p = netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
ggsave("Manifold_functional_cluster.pdf", p, width = 8, height = 6)
#netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

# 基于拓扑相似性的流行学习与分类
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
p = netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
ggsave("Manifold_structural_cluster.pdf", p, width = 8, height = 6)

## The end
saveRDS(cellchat, file = "cellchat.rds")

# 二、不同分组之间的配对分析####
# 配对分析必须保证细胞类型是一样的，才可以进行配对。如果 两个样本的细胞类型不一样又想进行配对分析时，可以用subset把两个样本的细胞类型取成一致的
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
rm(list = ls())
options(stringsAsFactors = FALSE)

## 创建cellchat对象
### 提取数据子集
scRNA <- readRDS("~/project/Integrate/scRNA.classified.rds")
scRNA$celltype <- scRNA$SingleR
table(scRNA$celltype)
Idents(scRNA) <- 'celltype'
scRNA <- subset(scRNA, idents = c('B cells','CD4+ T cells','CD8+ T cells','Dendritic cells','Monocytes','NK cells'))
scRNA$celltype <- as.factor(as.character(scRNA$celltype))
table(scRNA$orig.ident)
Idents(scRNA) <- 'orig.ident'
sco.til <- subset(scRNA, idents = c('HNC01TIL', 'HNC10TIL', 'HNC20TIL'))
sco.pbmc <- subset(scRNA, idents = c('HNC01PBMC', 'HNC10PBMC', 'HNC20PBMC'))

### 创建cellchat对象
cco.til <- createCellChat(sco.til@assays$RNA@data, meta = sco.til@meta.data, group.by = "celltype")
cco.pbmc <- createCellChat(sco.pbmc@assays$RNA@data, meta = sco.pbmc@meta.data, group.by = "celltype")
save(cco.til, cco.pbmc, file = "cco.rda")

# 2. 细胞通讯网络分析
# 2.1 数据准备和路径切换
dir.create("./Compare")
setwd("./Compare")
# load("../cco.rda")
# cco.pbmc <- setIdent(cco.pbmc, ident.use = "celltype")
# cco.til <- setIdent(cco.til, ident.use = "celltype")

# 2.2 分析样本cco.pbmc的细胞通讯网络
# cellchat不太稳定，identifyOverExpressedGenes常出错（不出现进度条就是出错了），重启Rstudio后再运算。
cellchat <- cco.pbmc
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cco.pbmc <- cellchat
saveRDS(cco.pbmc, "cco.pbmc.rds")

# 2.3 分析样本cco.til的细胞通讯网络
cellchat <- cco.til
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cco.til <- cellchat
saveRDS(cco.til, "cco.til.rds")

# 2.4 合并cellchat对象
cco.list <- list(pbmc=cco.pbmc, til=cco.til)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)

# 3. 可视化
# 3.1 所有细胞群总体观：通讯数量与强度对比
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("Overview_number_strength.pdf", p, width = 6, height = 4)

# 数量与强度差异网络图
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# save as Diff_number_strength_net.pdf

# 数量与强度差异热图
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
# save as Diff_number_strength_heatmap.pdf

# 细胞互作数量对比网络图
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}
# save as Counts_Compare_net.pdf

# 3.2 指定细胞互作数量对比网络图
par(mfrow = c(1,2))
s.cell <- c("CD4+ T cells", "CD8+ T cells", "Monocytes")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[2]))
# save as Counts_Compare_select.pdf 10*6.5

# 3.3 保守和特异性信号通路的识别与可视化
## 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strengh.pdf", p, width = 10, height = 6)

# 3.4 流行学习识别差异信号通路
# 这里function的图出不来，只有structural的图可以出来
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
#netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)
    
saveRDS(cellchat, "cellchat.rds")

# 3.5 细胞信号模式对比(总体信号模式对比)
library(ComplexHeatmap)
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_all.pdf  10*6

# 输出信号模式对比
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_outgoing.pdf  10*6

# 输入信号模式对比
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_incoming.pdf  10*6

# 3.6 特定信号通路的对比
# 网络图
pathways.show <- c("IL16") 
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(cco.list)[i]))
}
# save as Compare_IL16_net.pdf  10*6.5

# 热图
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# save as Compare_IL16_heatmap.pdf  12*6.5

# 和弦图
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "chord", pt.title = 3, title.space = 0.05,
                      vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(cco.list)[i]))
}
# save as Compare_IL16_chord.pdf  10*6.5

# 3.7 配体-受体对比分析
# 气泡图展示所有配体受体对的差异
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(1,2,3,6),  comparison = c(1, 2), angle.x = 45)
ggsave("Compare_LR_bubble.pdf", p, width = 12, height = 8)

# 气泡图展示上调或下调的配体受体对
p1 <- netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(1,2,3,6), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in TIL", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(1,2,3,6), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in TIL", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Compare_LR_regulated.pdf", pc, width = 12, height = 5.5)

# 和弦图
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]], sources.use = c(4,5), targets.use = c(1,2,3,6), signaling = "MHC-I", 
                       lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 20,
                       title.name = paste0("Signaling from Treg - ", names(cco.list)[i]))
}
# save as Compare_LR_chord.pdf  10*6.5
