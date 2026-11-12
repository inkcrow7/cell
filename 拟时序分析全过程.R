#install.packages(packageurl, repos=NULL, type ="source")
#install.packages(packageurl, repos=NULL, type ="source")
#install.packages("devtools")
#devtools::install_github("cole-trapnell-lab/monocle-release@develop")
#install.packages("monocle",version = "2.30.1")
# 创建CellDataSet
# 2.1将Seurat object中数据提取来创建
library(dplyr)
library(Seurat)
library(patchwork) # 用来做拼图的包，后面的p1|p1|p3在一张图上展示三个图就是这个包的功劳

# Load the PBMC dataset
pbmc.data<-Read10X(data.dir="D:/SG/单细胞/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw(non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") # 计算线粒体QC指标
# 提取细胞中基因数目大于200和小于2500的细胞，并且细胞中线粒体的比例小于5 （低--低质量细胞；高--可能是双细胞）
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method ="vst",nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features =all.genes)
pbmc <- RunPCA(pbmc,features =VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc,dims =1:10)
pbmc <- FindClusters(pbmc,resolution=0.5)
pbmc <- RunUMAP(pbmc,dims =1:10)

new.cluster.ids <- c("Naive CD4 T","Memory CD4 T","CD14+ Mono","B","CD8 T","FCGR3A+ Mono",
                    "NK","DC","Platelet")
names(new.cluster.ids)<- levels(pbmc)
pbmc<- RenameIdents(pbmc,new.cluster.ids)
pbmc[['cell_type']]<- pbmc@active.ident # 将注释结果添加到metadata
levels(pbmc)
saveRDS(pbmc,file="pbmc.rds") # 运行完上述代码，就会得到注释好的数据集

# monocle输入的是count矩阵(不建议使用data矩阵）
library(monocle) # igraph
pbmc<- readRDS("D:\\SG\\单细胞组学_9_9\\2.单细胞常规流程\\pbmc.rds") # 导入注释好的seurat对象(已注释)
# 提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix <- as(as.matrix(pbmc@assays$RNA$counts),'sparseMatrix')
# 提取表型信息到p_data(phenotype_data)里面
p_data<- pbmc@meta.data
#p_data$celltype<- pbmc@active.ident ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
# 提取基因信息 如生物类型、gc含量等
f_data<- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))
# expr_matrix的行数与f_data的行数相同(gene number)，expr_matrix的列数与p_data的行数相同(cell number)

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds<- newCellDataSet(expr_matrix,
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit = 0.5, # 原始计数（整数）：0.5；标准化数据（连续值）：0.1
                     expressionFamily = negbinomial.size()) # 输入的表达矩阵为UMI,一般适用于10x的数据；稀疏矩阵用negbinomial.size()
# tobit():适用于输入的表达矩阵为FPKM或者TPM；gaussianff():输入为log化后的FPKM或者TPM。
# FPKM/TPM值通常是对数正态分布的，而UMls或读计数使用负二项更好地建。要处理计数数据，需要将负二项分布指定为newCellDataSet的expressionFamily参数

# 方式二########################################################################
# 2.2直接读取表达矩阵来创建
library(data.table)
##读取数据
data <- fread("fpkm.txt",data.table=F,header =T)
pd <- fread("metadata.txt",data.table =F,header =T)
fd <- fread("gene annotations.txt",data.table =F,header = T)##创建
pd <- new("AnnotatedDataFrame", data =pd)
fd<- new("AnnotatedDataFrame", data =fd)
HSMM<- newCellDataSet(as.matrix(data),
                      phenoData =pd,featureData=fd,
                      expressionFamily = tobit())
# 如果数据量大，建议转化为稀疏矩阵
HSMM<-newCellDataSet(as(as.matrix(data),"sparseMatrix"),
                     phenoData =pd,featureData =fd,
                     expressionFamily = tobit())
####方式三####
importCDS(pbmc) # 如果我们想要从Seurat对象或SCESet中导入所有的插槽，我们可以设置参数'import_all'为TRUE。
#(默认为FALSE或只保留最小数据集)
#####

# 3估计size factor和离散值
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds) # 类似于seurat的数据归一化处理

# 4过滤低质量的细胞
cds <- detectGenes(cds, min_expr = 0.1) # 这一操作会在fData(cds)中添加-列num_cells_expressed
print(head(fData(cds))) #此时有13714个基因
expressed_genes <- row.names(subset(fData(cds),
      num_cells_expressed >= 10))#过滤掉在小于10个细胞中表达的基因，还剩11095个基因。

# 5细胞分类
# 建议先将细胞注释好再进行Monocle分析，不建议使用monocle做细胞分类。
# 无监督分析方法#####

# 6轨迹定义基因选择及可视化和构建轨迹
# 方法一：使用seurat选择的高变基因A
express_genes <- VariableFeatures(pbmc)
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)
# 方法二：使用clusters差异表达基因(一般选择这个方法)
deg.cluster <- FindAllMarkers(pbmc)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)
# 方法三：使用monocle选择的高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

# 方法四：dpFeature
# 这一步输入的expressed_genes来自于步骤4。后续分析使用的是该方法(半监督)
# 也可输入seurat筛选出的高变基因:expressed_genes<- VariableFeatures(pbmc) 
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~cell_type",cores=1)#~后面是表示对谁做差异分析的变量，理论上可以为pdata的任意列名
# 这步较慢
head(diff)
# 差异表达基因作为轨迹构建的基因,差异基因的选择标准是qva1<0.01,decreasing=F表示按数值增加排序
deg <- subset(diff,qval < 0.01)#选出2829个基因,或者（num_cells_expressed >= 10）
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

##差异基因的结果文件保存
write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)

##轨迹构建基因可视化
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
# 这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setorderingfilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
# setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
# 出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)

# gene数太多的话也可以选择top基因
#ordergene <- row.names(deg)[order(deg$qval)][1:400]

# 6.2降维
cds <- reduceDimension(cds,max_components=2, # 这步较慢
                       reduction_method  = 'DDRTree')
# 6.3拟时间轴轨迹构建和在拟时间内排列细胞
cds <- orderCells(cds)
# 使用root_state参数可以设置拟时间轴的根,如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds,root_state=5)#把state5设成拟时间轴的起始点

# 可视化：根据cds@phenoData@data中的表型信息(metadata)对细胞上色
pdf("train.monocle.pseudotime.pdf",width=7,height=7)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE)
dev.off()

# 以细胞类型上色
pdf("train.monocle.celltype.pdf",width=7,height =7)
plot_cell_trajectory(cds,color_by="cell_type", size=1,show_backbone=TRUE)
dev.off()

# 以细胞状态上色
pdf("train.monocle.state.pdf",width=7,height=7)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
dev.off()

# 按照seurat分群排序细胞
pdf("seurat.clusters.pdf",width=7,height =7)
plot_cell_trajectory(cds, color_by = "seurat_clusters")
dev.off()

# 以细胞状态上色（拆分）“分面”轨迹图，以便更容易地查看每个状态的位置。
pdf("train.monocle.state.faceted.pdf",width= 10,height = 7)
plot_cell_trajectory(cds, color_by = "State")+ facet_wrap("~State",nrow = 1)
dev.off()

library(ggsci)
p1=plot_cell_trajectory(cds, color_by="cell_type")+ scale_color_npg()
p2=plot_cell_trajectory(cds, color_by="State")+ scale_color_nejm()
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#7f7f7f", "#bcbd22")
p3=plot_cell_trajectory(cds, color_by = "State")+ scale_color_manual(values = colour)
p1|p2|p3

# 觉得这种轨迹不太直观的也可以画成树形图
p1<-plot_cell_trajectory(cds,x=1,y=2,color_by="celltype")+
  theme(legend.position='none',panel.border = element_blank())+ #去掉第一个的legend
  scale_color_manual(values = colour)
p2<-plot_complex_cell_trajectory(cds,x=1,y=2,color_by="celltype")+
  scale_color_manual(values = colour)+
  theme(legend.title = element_blank())
p1|p2

# 还可以画沿时间轴的细胞密度图
library(ggpubr)
df <- pData(cds)
## pData(cds)取出的是cds对象中cds@phenoData@data的内容
View(df)
ggplot(df, aes(Pseudotime, colour = cell_type, fill=cell_type))+
  geom_density(bw=0.5,size=1,alpha=0.5)+
  theme_classic2()

#手动设置颜色
ClusterName_color_panel <- c(
  "Naive CD4 T" = "#DC143C","Memory CD4 T" = "#0000FF","CD14+ Mono" = "#20B2AA",
  "B" = "#FFA500","CD8 T" = "#9370DB","FCGR3A+ Mono" = "#98FB98",
  "NK" = "#F08080","DC" = "#0000FF","Platelet" = "#20B2AA"
)
ggplot(df, aes(Pseudotime,colour = cell_type,fill=cell_type)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+
  theme_classic2()+scale_fill_manual(name = "",
                                     values = ClusterName_color_panel)+
  scale_color_manual(name = "",values = ClusterName_color_panel)

# 提取感兴趣的细胞（进行后续分析）
#比如对state7的细胞感兴趣
pdata <- Biobase::pData(cds)
s.cells <- subset(pdata, State=="3")%>% rownames()
save(s.cells,file ="Monocle state3.rda")

write.csv(pData(cds),"pseudotime.csv")
save(cds,file ="cds.rda")

# 指定基因的可视化
##选择前4个top基因并将其对象取出
keygenes <- head(ordergene,4)
cds_subset <- cds[keygenes,]
##可视化:以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
ggsave("Genes_pseudotimeplot.pdf", plot= plotc, width= 16, height = 8)

# 指定基因
s.genes <- c("SELL","CCR7","IL7R","CD84","CCL5","S100A4")
p1 <- plot_genes_jitter(cds[s.genes,],grouping = "State", color_by = "State")
p2 <- plot_genes_violin(cds[s.genes,],grouping = "State", color_by = "State")
p3<- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
plotc<- p1|p2|p3
ggsave("Genes_Jitterplot.pdf",plot = plotc, width= 16,height = 8)

# 拟时序展示单个基因表达量
colnames(pData(cds))
pData(cds)$CCL5 =log2(exprs(cds)['CCL5',]+1)
p1=plot_cell_trajectory(cds, color_by="CCL5") + scale_color_gsea()
pData(cds)$S100A4=log2(exprs(cds)['S100A4',]+1)
p2=plot_cell_trajectory(cds, color_by="S100A4")+ scale_color_gsea()
library(patchwork)
p1+p2

#8. 寻找拟时相关的基因（拟时差异基因）
#使用回归算法(注意：不要使用多核运算，经常会出现警告)

# 寻找拟时差异基因（qvalue体现基因与拟时的密切程度）绘制热图
#这里是把排序基因(ordergene)提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1,
                                  fullModelFormulaStr ="~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)]#把gene放前面，也可以不改
write.csv(Time_diff,"Time_diff_all.csv",row.names =F)
Time_genes <- Time_diff %>% pull(gene_short_name)%>% as.character()
p=plot_pseudotime_heatmap(cds[head(Time_genes,650),], 
                          num_clusters=4, 
                          show_rownames=T, 
                          return_heatmap=T,
                          #hmcols = colorRampPalette(c("navy","white","firebrick3"))(100)
  )

ggsave("Time_heatmapAll.pdf",p,width=5,height = 10)

# 把每个cluster的基因单独提出来做分析
p$tree_row
# Ca11:
# hclust(d=d, method = method)
# Cluster method :ward.D2
# Number of objects:2829
clusters <- cutree(p$tree_row,k=4) # 提取第四个cluster的基因
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_clusters"
table(clustering)
# 1    2    3    4 
# 570 1031  506  722 
write.csv(clustering,"Time_clustering_all.csv",row.names =F)

# 拟时差异基因热图绘制（提取了前100个）
Time_genes <- top_n(Time_diff,n= 100,desc(qval))%>% pull(gene_short_name) %>% as.character()
p= plot_pseudotime_heatmap(cds[Time_genes,], 
                           num_clusters=4,
                           show_rownames=T,
                           return_heatmap=T)
ggsave("Time_heatmapTop100.pdf",p,width=5,height = 10)

# 显著差异基因按热图结果排序并保存
hp.genes<- p$tree_row$labels[p$tree_row$order]
Time_diff_sig<- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(Time_diff_sig, "Time_diff_sig.csv",row.names = F)

# 手动选择基因来绘制热图，查看其表达模式
marker_genes <- row.names(subset(fData(cds),
                               gene_short_name %in% c("MEF2C","MEF2D","'MYF5",
                                                     "ANPEP","PDGFRA" ,"MYOG",
                                                     "TPM1","TPM2","MYH2",
                                                     "MYH3","NCAM1","TNNT1",
                                                     "CDK1","TNNT2" ,"TNNC1",
                                                     "CDK2","CCNB1","CCNB2",
                                                     "CCND1","CCNA1","ID1")))
diff_test_res <- differentialGeneTest(cds[marker_genes,],
                                    fullModelFormulaStr ="~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 6,
                        cores = 1,
                        show_rownames =T)

# 单细胞轨迹的“分支”分析
plot_cell_trajectory(cds, color_by ="State")

# BEAM进行统计分析将ordergene换为[head(Time_genes,650),]
BEAM_res <- BEAM(cds[head(Time_genes,650),],branch_point = 1,cores = 2,progenitor_method = "duplicate")
#这里用的是ordergene，也就是第六步dpFeature找出来的基因。如果前面用的是seurat的mark基因，记得改成express_genes
#BEAM_res<- BEAM(cds,branch point=1,cores=2)#对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name","pval","qval")]
head(BEAM_res)
#           gene_short_name         pval         qval
# CD79A               CD79A 2.757676e-73 7.782161e-70
# TCL1A               TCL1A 1.574889e-65 2.222168e-62 
# IGLL5               IGLL5 2.356778e-64 2.216942e-61
# S100A9             S100A9 1.504319e-58 1.061297e-55
# S100A8             S100A8 6.028175e-57 3.402302e-54
# LINC00926       LINC00926 3.180527e-55 1.495908e-52
write.csv(BEAM_res,"BEAM_res.csv", row.names = F)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-4)),],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T) #有632个gene，太多了

# 选前100个基因可视化
BEAM_genes <- top_n(BEAM_res,n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p<- plot_genes_branched_heatmap(cds[BEAM_genes,], branch_point = 1,
                                num_clusters =3,show_rownames =T,return_heatmap =T)
ggsave("BEAM_heatmap.pdf",p$ph_res,width= 6.5,height = 10)

#显著差异基因(top100)按热图结果排序并保存
##如果要所有的差异基因，就把前面所632个基因的热图存为p
hp.genes<-p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig<- BEAM_res[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(BEAM_sig,"BEAM_sig.csv",row.names =F)

# 选择上面热图中或显著差异基因中感兴趣的基因进行可视化
head(BEAM_sig)
genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c( "LST1", "CEBPD", "TYROBP")))
plot_genes_branched_pseudotime(cds[genes,],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1)
save.image(file = "拟时序分析.RData")
