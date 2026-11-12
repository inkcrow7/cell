#1.在电脑说安装JAGS（前提条件），根据链接https://sourceforge.net/projects/mcmc-jags/下载JAGS 并安装(和R的安装路径一致)
#2.BiocManager::install("rjags")
#BiocManager::install("infercnv")
library(infercnv) 
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
                                    annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", 
                                    package = "infercnv"),
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=tempfile(),
                             cluster_by_groups=TRUE,
                             denoise=TRUE,
                             HMM=TRUE)


# 实例
rm(list = ls())
library(infercnv)
library(Seurat)
library(gplots)
library(ggplot2)
library(qs)
sce <- qread("sce.all.qs")

#检查一下自己导入进来的数据
DimPlot(sce,reduction = 'umap',
        label = TRUE,pt.size = 0.5) +NoLegend() # pt.size调整点的大小

dir.create("11-infercnv")
setwd("./11-infercnv")

table(Idents(sce))
table(sce@meta.data$seurat_clusters)
table(sce@meta.data$orig.ident)

# 提取部分细胞
sce@meta.data
sce <- sce[,sce$celltype %in% c("plasma","fibroblasts","epithelial/cancer cells")]

# 文件制作1：表达量矩阵
# 除了用GetAssayData函数其实也可以直接sce@assays$RNA@count即可
dat <- GetAssayData(sce,layer = 'counts',assay = 'RNA')
dat[1:4,1:4]

# 文件制作2：样本的描述
groupinfo <- data.frame(v1 = colnames(dat),v2 = sce@meta.data$celltype)
head(groupinfo)
#                                  v1                      v2
# 1 GSM5688706_WGC_AACACACGTGCTTCAA-1 epithelial/cancer cells
# 2 GSM5688706_WGC_AACACACTCCCGTTCA-1 epithelial/cancer cells
# 3 GSM5688706_WGC_AACACACTCTGCTAGA-1                  plasma
# 4 GSM5688706_WGC_AACCATGAGGCACAAC-1                  plasma
# 5 GSM5688706_WGC_AACCATGCAAAGCTCT-1 epithelial/cancer cells
# 6 GSM5688706_WGC_AACCCAAGTGCCGTAC-1                  plasma

# 文件制作3：基因在染色体中的坐标
library(AnnoProbe)
#annoGene 函数返回一个数据框，包含输入基因的详细注释信息。
#注释信息可能包括基因名称、染色体位置、基因描述等。
geneInfor <- annoGene(rownames(dat),"SYMBOL","human") #物种需要小心哦
# 使用逻辑条件来移除包含 chrM, chrX, 和 chrY 的行
geneInfor <- geneInfor[!geneInfor$chr %in% c("chrM", "chrX", "chrY"), ]
#提取chr后后面的数字并转化为num，从而按这个num排序
#sub函数用于把chr替换为空
geneInfor$chr_num <- as.numeric(sub("chr", "", geneInfor$chr))
colnames(geneInfor)

#with函数的作用是简化写法，可问gpt
geneInfor <- geneInfor[with(geneInfor,order(chr_num,start)),c(1,4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]

length(unique(geneInfor[,1]))
head(geneInfor)
#        SYMBOL  chr  start    end
# 23 AP006222.2 chr1 266855 268655
# 51     FAM87B chr1 817371 819837
# 53  LINC00115 chr1 826206 827522
# 55     FAM41C chr1 868071 876903
# 62     SAMD11 chr1 923928 944581
# 63      NOC2L chr1 944203 959309

#保留行名在geneInfor第一列中存在的行。
dat <- dat[rownames(dat) %in% geneInfor[,1],] 
#match(x, table)：match函数返回x中每个元素在table中的位置索引。
#获得位置后对dat进行重新排序使其跟geneInfor中的顺序一致
dat <- dat[match(geneInfor[,1],rownames(dat)),]
dim(dat)
head(groupinfo)
dat[1:4,1:4]
table(groupinfo$v2)
dim(groupinfo)

#为了节约计算机资源，直接抽样
#如果真实项目请不要用抽样
kp <- sample(1:nrow(groupinfo),500)
groupinfo <- groupinfo[kp,]
dat <- dat[,kp]

# 保存矩阵文件
#expFile:是一个变量，存储写入的文件的文件名或路径。在这里文件名是expFile.txt。
expFile <- 'expFile.txt'#定义输出文件名
colnames(dat) <- gsub("-", "_", colnames(dat))
write.table(dat, file = expFile, sep = '\t', quote = F)
# 保存细胞注释文件
groupFiles <- 'groupFiles.txt'
groupinfo$v1 <- gsub("-", "_", groupinfo$v1)
write.table(groupinfo,file = groupFiles, sep = '\t',
            quote = F, col.names = F, row.names = F)
# 保存基因排序文件
head(geneInfor)
geneFile <- 'geneFile.txt'
write.table(geneInfor, file = geneFile, sep = '\t',
            quote = F, col.names = F, row.names = F)

# ref_group_names选择已知的非肿瘤细胞即可，我们这里选择了plasma细胞；
# 如果没有参考细胞，可以将 ref_group_names 设置为 NULL，在这种情况下，将使用所有细胞的平均信号来定义基线。
#创建对象，请注意文件需要一一对应哦！
#ref_group_names 这里的细胞是正常对照，然后跟其他的细胞比较
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = expFile,
                                     annotations_file = groupFiles,
                                     delim = "\t",
                                     gene_order_file = geneFile,
                                     ref_group_names = c("plasma","fibroblasts") #这里用plasma和fibroblasts做参照
)

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1, #smart-seq适应于1,10X适应于0.1
                              out_dir = "infercnv_output", # dir is auto
                              # 是否根据细胞注释文件的分组
                              # 对肿瘤细胞进行分组
                              # 影响read.dendrogram, 如果有多个细胞类型，且设置为TRUE，
                              # 后续的read.dendrogram无法执行
                              cluster_by_groups = TRUE, #是否根据细胞类型(由细胞注释文件中定义)
                              hclust_method = "ward.D2",# ward.D2 方法进行层次聚类
                              analysis_mode = "subclusters", # 默认是samples，推荐是subclusters
                              denoise = TRUE, # 去噪音
                              HMM = F,  ##特别耗时间,是否要去背景噪音
                              plot_steps = F, #不在每个步骤后生成图形。
                              leiden_resolution = "auto", #可以手动调参数
                              #num_threads = 16# 16线程工作，加快速度，一般电脑不建议哈
)

