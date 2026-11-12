# install.packages('Seurat')
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/home/taojj/单细胞组学_9_9/2.单细胞常规流程/filtered_gene_bc_matrices/hg19")
# 使用原始数据（非标准化数据）初始化Seurat对象
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# min.cells = 3表示一个基因至少在3个细胞中表达才会被保留;min.features = 200表示一个细胞至少检测到200个基因才会被保留
pbmc

# 使用函数计算线粒体QC指标，该函数可计算来自一组特征的计数百分比PercentageFeatureSet()，正则表达式匹配以MT-开头的基因名
# QC主要的指标有nCount_RNA(测到的转录本分子总数（UMI数目）测到的基因总数)和nFeature_RNA（每个细胞中检测到的基因数量）以及"percent_mito"（表示细胞中线粒体基因的比例）这三个指标。
# 可选"percent_ribo"（核糖体基因比例）和"percent_hb"（红血细胞基因比例）两个指标。
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 可视化
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter 通常用于可视化特征与特征之间的关系，但也可用于对象计算的任何内容，如对象元数据中的列，PC 分数等。
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# nFeature_RNA > 200是为了去除低质量的细胞；nFeature_RNA < 2500则是过滤掉高基因数的细胞；percent.mt < 5则是控制线粒体基因比例，过高的比例通常指示细胞凋亡或损伤。
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) # percent.mt严格小于5，宽松一点小于10。特殊情况下（肝脏）小于20

# 标准化数据
# 采用全局缩放归一化方法“LogNormalize”，该方法通过总表达式对每个单元格的特征表达式测量值进行归一化，
# 将其乘以比例因子（默认为 10,000），然后对结果进行对数转换。在Seurat v5 中，标准化值存储在pbmc[["RNA"]]$data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc <- NormalizeData(pbmc)

# 识别高度可变的特征（特征选择）
# 默认情况下，我们为每个数据集选取 2000 个特征。这些将用于下游分析，如 PCA
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 确定高表达的前十个基因
top10 <- head(VariableFeatures(pbmc), 10)

# 绘制带标签和不带标签的变量特征图
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.new(width=12, height=8, unit="in")  # 创建新窗口（12英寸宽 x 8英寸高）
plot1 + plot2

# 缩放数据（调整基因表达值，使得不同细胞之间的表达水平具有可比性，同时为后续的降维（如PCA）做准备）
# 应用线性变换（“缩放”），这是 PCA 等降维技术之前的标准预处理步骤
all.genes <- rownames(pbmc)
# 改变每个基因的表达，使细胞间的平均表达量为 0 ；缩放每个基因的表达，使细胞间的方差为 1
pbmc <- ScaleData(pbmc, features = all.genes) # 结果存储在pbmc[[“RNA”]]@scale.data

# 主成分分析（PCA）
# 线性降维
# Seurat 提供了几种有用的方法来可视化定义 PCA 的单元格和特征，包括VizDimReduction()、DimPlot()和DimHeatmap()
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) # 选择第一主成分，正负各250个细胞

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) # # 选择1:15主成分，正负各250个细胞

# 确定数据集的“维度” 
pbmc<-JackStraw(pbmc,num.replicate = 100)
pbmc<-ScoreJackStraw(pbmc,dims = 1:20)

# 可视化处理
JackStrawPlot(pbmc,dims = 1:15) # 选择P值显著的主成分，13到14变为不显著，因此不选

# 根据每个成分（函数）解释的方差百分比对主成分进行排名
ElbowPlot(pbmc)

# 聚类细胞
# 建立KNN图，并基于其局部领域中的共享重叠细化任意两个单元之间的边权重
pbmc<-FindNeighbors(pbmc,dims = 1:10)
# 对细胞进行聚类，应用模块化优化技术，Louvain算法（默认）或SLM
# 以迭代方式将细胞分组在一起，目标是优化标准模块化函数
pbmc<-FindClusters(pbmc,resolution = 0.5) # resolution参数控制聚类的粒度（granularity）。
#resolution值越大，得到的聚类数越多（更细的粒度）；值越小，聚类数越少（更粗的粒度）
head(Idents(pbmc),5)

# 非线性降维（UMAP/tSNE）
pbmc<-RunUMAP(pbmc,dims = 1:10)
DimPlot(pbmc,reduction = "umap") # pt.size = 2.5
# DimPlot(kid, reduction = "umap", group.by = "ident", pt.size=0.5)
# DimPlot(kid, reduction = "umap", group.by = "orig.ident", pt.size=0.5,split.by = 'orig.ident')

# tSNE
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")
# 运行t-SNE（带参数优化）
pbmc <- RunTSNE(
  object = pbmc,
  dims = 1:10,
  perplexity = 30,          # 根据细胞数量调整
  max_iter = 1500,          # 增加迭代保证收敛
  seed.use = 42,            # 固定随机种子
  check_duplicates = FALSE  # 无重复细胞ID时关闭
)
# 寻找差异表达的特征（簇生物标志物）
# findmarkers为所有集群自动执行此过程，也可以测试集群组之间的对比，或针对所有单元格进行测试
# 默认情况下，ident.1与所有其他细胞相比，他识别单个簇的阳性和阴性标记。
# min.pct参数要求在两组细胞中的任何一组中以最小百分比检测到一个特征
# 而thresh.test参数要求一个特征在两组之间差异表达（平均）一定量
# 寻找cluster2的所有markers
cluster2.markers<-FindMarkers(pbmc,ident.1 = 2,min.pct = 0.25)
head(cluster2.markers,n=5)

# 寻找cluster5中与cluster0和cluster3n不同的所有markers
cluster5.markers<-FindMarkers(pbmc,ident.1 = 5,ident.2=c(0,3),min.pct = 0.25)
head(cluster5.markers,n=5)

# 找出每个细胞簇的标记物，与所有剩余的细胞进行比较，只报告阳性细胞 
pbmc.markers<-FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 0.25,logfc.threshold =0.25 ) #only.pos = TRUE意味着只保留在目标cluster中上调的基因，排除下调的。
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n=2,order_by = avg_log2FC) %>% 
  print(n = 18)

# 每个亚群topn特征基因
topn_genelist <- pbmc.markers %>%
     filter(pct.1 >= 0.5, cluster != "Unknown") %>%
     group_by(cluster) %>% 
     slice_max(n = 5, order_by = avg_log2FC)

# Dotplot展示
library(ggplot2)
p <- DotPlot(object = pbmc, features = unique(topn_genelist$gene), 
             cols = c("#C0C0C0", "#E41C12")) +
     scale_color_gradient2(low = "blue", high = "red", space = "Lab", mid = "#C0C0C0") +
     viridis::scale_color_viridis() +
     theme(
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
         axis.text = element_text(color = "black")
    )

# 还有用于可视化标记表达的工具，Vlnplot(显示跨集群的表达概率分布)和FeaturePlot()(在tSNE或PCA图上可视化特征表达)是最常用的可视化，建议探索RidgePlot(),CellScatter（），和DotPlot()作为查看数据集的其他方法
VlnPlot(pbmc,features = c("MS4A1","CD79A"))
# 可以绘制行数据
VlnPlot(pbmc,features = c("NKG7","PF4"),slot = "counts",log = TRUE)

FeaturePlot(pbmc,features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                              "CD8A"))

#DoHeatmap()为给定的细胞和特征生成一个表达热图。在这种情况下，我们为每个集群绘制前 20 个标记（或所有标记，如果小于 20）。
pbmc.markers%>%
  group_by(cluster)%>%
  top_n(n=10,wt=avg_log2FC)->top10
DoHeatmap(pbmc,features = top10$gene)+NoLegend()

# 将细胞类型标识分配给集群
new.cluster.ids<-c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                   "NK", "DC", "Platelet")
names(new.cluster.ids)<-levels(pbmc)
pbmc<-RenameIdents(pbmc,new.cluster.ids)
DimPlot(pbmc,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()
saveRDS(pbmc,file = "pbmc3k_final.rds")

# 提取某一细胞类型
table(pbmc@active.ident)
sce1 = pbmc[, pbmc@active.ident %in% c( 'Naive CD4 T' )]
