#rm(list=ls())
library(Seurat)
library(gplots)
library(ggplot2)
#load('sce.markers.all_10_celltype.Rdata')

library(clusterProfiler)
library(org.Hs.eg.db)
ids=bitr(pbmc.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
pbmc.markers=merge(pbmc.markers,ids,by.x='gene',by.y='SYMBOL')
View(sce.markers)

# 富集分析绘图
GO <- enrichGO(  
  markers$ENTREZID,  
  OrgDb = 'org.Hs.eg.db',  
  keyType = "ENTREZID",  
  ont = "ALL",  
  pvalueCutoff = 0.01,  
  pAdjustMethod = "BH",  
  qvalueCutoff = 0.05,  
  minGSSize = 50,  
  maxGSSize = 500,  
  readable = TRUE)

KEGG <- enrichKEGG(  
  markers$ENTREZID, 
  organism = "hsa",   
  keyType = "kegg",  
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10, 
  maxGSSize = 500, 
  qvalueCutoff = 0.05,  
  use_internal_data = FALSE) 
##GO富集分析绘图
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free") +  
  scale_y_discrete(labels = function(x)     
    stringr::str_wrap(x, width = 80)) ##避免字体重合
#dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")

##KEGG富集分析绘图
barplot(KEGG,showCategory = 15,title = 'KEGG Pathway')

##可视化信号通路
library("pathview")
gene$ENTREZID <- as.numeric(gene$ENTREZID )
pathview(gene.data = markers$ENTREZID,          
         pathway.id = "hsa05171",          
         species = "hsa",         
         out.suffix = "hsa05171.2layer",         
         kegg.native = T, same.layer = F)

## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(pbmc.markers$ENTREZID, pbmc.markers$cluster) 
## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(xx,label_format = 100)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

## GO
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)
p <- dotplot(xx,label_format = 100)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

# 差异分析后的GO以及KEGG分析
# 读取数据
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(remotes)
# BiocManager::install("MAST")
library(MAST)
load(file = 'basic.sce.pbmc.Rdata')
sce=pbmc
sce = sce[, Idents(sce) %in% c("FCGR3A+ Mono", "CD14+ Mono")] # 挑选细胞

# 差异分析
deg=FindMarkers(object = sce,
                ident.1 = 'FCGR3A+ Mono',
                ident.2 = 'CD14+ Mono', 
                test.use='MAST')  ## MAST在单细胞领域较为常用
head(deg)
save(deg,file = 'deg-by-MAST-for-mono-2-cluster.Rdata')

# 火山图
# 单细胞的logFC不像普通测序那样大于1很多，但是p值小于0.05贼多，所以火山图参考一篇CNS 文章，只画了p值的分界线
degdf <- deg
degdf$symbol <- rownames(deg)
logFC_t=0
P.Value_t = 1e-28
degdf$change = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < 0,"down",
                      ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 0,"up","stable"))
ggplot(degdf, aes(avg_log2FC,  -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=3.5, aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()

# 功能注释
## 获取上下调基因
gene_up=rownames(deg[deg$avg_log2FC > 0,])
gene_down=rownames(deg[deg$avg_log2FC < 0,])
## 把SYMBOL改为ENTREZID
library(org.Hs.eg.db)
gene_up=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                   keys = gene_up,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))
gene_down=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                     keys = gene_down,
                                                     columns = 'ENTREZID',
                                                     keytype = 'SYMBOL')[,2]))
library(clusterProfiler)
## 以上调基因为例，下调基因同理
## KEGG
gene_up <- unique(gene_up)
kk.up <- enrichKEGG(gene = gene_up,
                    organism = "hsa",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
dotplot(kk.up)
barplot(kk.up)
# 结果存储
kk.up <- DOSE::setReadable(kk.up,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID") # 将ENTREZID转换为genesymbol
write.csv(data.frame(ID=rownames(kk.up@result),kk.up@result),
          file = "KEGG_enrich_result_kk.up.csv")

## GO
go.up <- enrichGO(gene = gene_up,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP" ,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.99,
                  qvalueCutoff = 0.99,
                  readabl = TRUE)
dotplot(go.up)
barplot(go.up)
write.csv(data.frame(ID=row.names(go.up@result),go.up@result),
          file = "GO_result.csv")

# 差异分析后的GSEA分析
## 上一步差异分析得到差异基因列表deg后取出，p值和log2FC
nrDEG = deg[,c('avg_log2FC', 'p_val')]
colnames(nrDEG)=c('log2FoldChange','pvalue') ##更改列名
library(org.Hs.eg.db)
library(clusterProfiler)
## 把SYMBOL转换为ENTREZID，可能有部分丢失
gene <- bitr(rownames(nrDEG),     
             fromType = "SYMBOL",     
             toType =  "ENTREZID",    
             OrgDb = org.Hs.eg.db)
## 基因名、ENTREZID、logFC一一对应起来
gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))]
## 构建genelist
geneList=gene$logFC
names(geneList)=gene$ENTREZID 
geneList=sort(geneList,decreasing = T) # 降序，按照logFC的值来排序
## GSEA分析
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 10,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
kk_gse=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
sortkk <- kk_gse[order(kk_gse$enrichmentScore, decreasing = T),]
library(enrichplot)
gseaplot2(kk_gse, 
          "hsa04510", 
          color = "firebrick",
          rel_heights=c(1, .2, .6))

## 展示排名前四的通路
gseaplot2(kk_gse, row.names(sortkk)[1:4])
