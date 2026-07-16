
# 安装R包
# install.packages("ggplot2")
# install.packages("patchwork")

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
setwd("F:/4scRNA/1code/4QC")

#一定要用fread，否则其他函数要很久
pbmc <-fread("GSE141445/data.raw.matrix.txt", sep="\t")
pbmc[1:5,1:5]

#直接转化为行名
pbmc <-  column_to_rownames(pbmc,"V1")
pbmc[1:5,1:5]

#创建seurat对象
pbmc <- CreateSeuratObject(pbmc, min.features = 200, min.cells = 3)

#将Layers融合
pbmc = JoinLayers(pbmc)

#查看
table(pbmc@meta.data$orig.ident)

table(str_split(colnames(pbmc),'-',simplify = T)[,2])
#增加或修改meta.dada信息
pbmc <- AddMetaData(object = pbmc, 
                    metadata = str_split(colnames(pbmc),'-',simplify = T)[,2],   
                    col.name = "orig.ident") 

table(pbmc@meta.data$orig.ident)

#%in% 判断前面一个向量内的元素是否在后面一个向量中，返回布尔值。
table(pbmc@meta.data$orig.ident %in% c("1","2","3"))

dim(pbmc)

#提取
pbmc = subset(pbmc,orig.ident %in% c("1","2","3"))

dim(pbmc)

#查看
table(pbmc@meta.data$orig.ident)

#质量控制
#线粒体基因比例
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(pbmc))
pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc, features=HB.genes) 

#查看相关性
FeatureScatter(pbmc, "nCount_RNA", "percent.mt", group.by = "orig.ident")
FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")

#查看质控指标
#设置绘图元素
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group = "orig.ident"
#质控前小提琴图
plots = list()
for(i in c(1:length(plot.featrures))){
  plots[[i]] = VlnPlot(pbmc, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)  
violin
#保存
ggsave("1vlnplot_before_qc.pdf", plot = violin, width = 14, height = 8) 
dim(pbmc)

#设置质控指标
quantile(pbmc$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(pbmc$nFeature_RNA, seq(0.9, 1, 0.01))
#plots[[1]] + geom_hline(yintercept = 500) + geom_hline(yintercept = 4500)
quantile(pbmc$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(pbmc$nCount_RNA, seq(0.9, 1, 0.01))
#plots[[2]] + geom_hline(yintercept = 22000)
quantile(pbmc$percent.mt, seq(0.9, 1, 0.01))
#plots[[3]] + geom_hline(yintercept = 20)
quantile(pbmc$percent.HB, seq(0.9, 1, 0.01))
#plots[[4]] + geom_hline(yintercept = 1)

#设置质控标准
#基因
minGene=300
maxGene=10000
#counts
minUMI=600
#线粒体
pctMT=10
#血细胞
pctHB=1

#数据质控并绘制小提琴图
pbmc <- subset(pbmc, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                 nCount_RNA > minUMI & percent.mt < pctMT & percent.HB < pctHB)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(pbmc, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
violin

#保存
ggsave("2vlnplot_after_qc.pdf", plot = violin, width = 14, height = 8) 
dim(pbmc)

#保存
saveRDS(pbmc,"1pbmc_qc.rds")

#读取
pbmc = readRDS("1pbmc_qc.rds")
