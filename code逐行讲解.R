# 安装加载所需要的R包（BiocManager）
#if (!requireNamespace("BiocManager", quietly = TRUE))  # 本行已注释，不执行；requireNamespace() 可检查是否已安装该包，quietly=TRUE 减少提示。当前下一行是直接安装，并不受此条件控制。
install.packages("BiocManager")  # 从 CRAN 安装 BiocManager，用于管理 Bioconductor 软件包；当前这行已启用，通常首次安装时运行即可。
library(GEOquery)  # 加载 GEOquery 包，用于读取 GEO 的系列矩阵文件（series matrix）。
BiocManager::install("GEOquery")  # 安装或更新 Bioconductor 的 GEOquery 包；当前这行已启用，且位于 library(GEOquery) 之后，首次使用应先确保包已安装。
library(stringr)  # 加载字符串处理包，提供 str_detect()、str_replace() 等函数（本脚本后续主要用基础函数，此包属预加载）。
library(ggplot2)  # 加载 ggplot2，本脚本所有正式图形（PCA、火山图）都基于它。
library(dplyr)  # 加载 dplyr，提供 %>% 管道、mutate()、case_when() 等数据框操作函数，用于整理分组信息。
library(ggpubr)  # 加载 ggpubr，用于发表级图形封装与组间显著性标注（本脚本主要作为工具备用）。
library(grid)  # 加载 grid 底层图形系统，提供 unit()、gpar() 等，热图注释和标签尺寸设置要用到。
library(ggrepel)  # 加载 ggrepel，提供 geom_text_repel()，使散点标签自动避让、不互相重叠。
library(limma)  # 加载 limma 包，本脚本差异表达分析的核心（线性模型 + 经验贝叶斯）。

# 数据下载与表达矩阵提取--------------------------------------------------------


# 设置工作路径
setwd("~/新转录组基础班20260724/文章复现2026724/2. 数据下载")  # 把工作目录切换到分析文件夹；~ 代表用户主目录。之后所有相对路径的读写都以该目录为基准。换机器运行时需要改成自己的路径。
getwd()  # 打印当前工作目录，用于确认 setwd() 是否生效。

# 数据下载与提取
library(GEOquery)  # 再次加载 GEOquery（重复加载无副作用，便于单独复制本段运行）。
eSet = getGEO("GSE179156",destdir = '.',getGPL = F, AnnotGPL = F)  # 下载并解析 GSE179156 的 Series Matrix，通常返回包含 ExpressionSet 的列表；destdir="." 保存到工作目录，getGPL=F 不另外获取平台记录，AnnotGPL=F 不请求专用 AnnotGPL 注释。

exp = exprs(eSet[[1]])  # 从列表的第一个 ExpressionSet 提取表达矩阵，行为探针、列为样本；若研究包含多个平台，应先确认选择的是目标平台。
exp = as.data.frame(exp)  # 把矩阵转成数据框，方便后续用数据框的方式改列、合并。
cli = pData(eSet[[1]])  # 提取样本表型信息（phenotype data），即每个样本的标题、分组、平台等字段，一行一个样本。

exp[exp <= 0] <- NaN  # 把不大于零的数值替换为 NaN，避免取对数时产生负无穷或无效值；这一步不代表所有缺失值都由非正数造成。
exp1 <- log2(exp)  # 对数值逐个取 log2，压缩表达范围；应先确认输入尚未取对数，不能对已是 log2 尺度的数据机械重复转换。
exp1 <- na.omit(exp1)  # 删除任意样本列中含 NA 或 NaN 的探针行；是按整行剔除，不是只删除单个缺失元素。

# 基因探针和symbol转换
# BiocManager::install("hgu133plus2.db")  # 网页中搜索查看GPL对应的Bioconductor注释包  # 被注释掉的安装语句：注释包按平台选择，需在 Bioconductor 网站搜索 GPL 编号对应的 .db 包；本数据平台对应 hgu133plus2.db。
library(hgu133plus2.db)  # 加载人全基因组 U133 Plus 2.0 芯片的注释数据库包。
probe2symbol = toTable(get("hgu133plus2SYMBOL"))  # get("hgu133plus2SYMBOL") 取出「探针 → symbol」的映射对象，toTable() 把它转成两列数据框（probe_id、symbol）。
rownames(probe2symbol) = probe2symbol$probe_id  # 把 probe_id 列设为行名，为下一步按行名合并做准备。

exp2 = merge(probe2symbol,exp1,by="row.names")  # 按行名中的探针 ID 合并注释表和表达矩阵；merge() 默认是内连接，只保留两边都有的探针，并新增 Row.names 列。
exp3 = aggregate(x=exp2, by=list(exp2$symbol), FUN=mean, na.rm=T)  # 按 symbol 分组，对同一基因的多个探针表达求均值；na.rm=T 忽略缺失。当前传入整张 exp2，字符型注释列也可能触发求均值警告，随后会删除这些列。
rownames(exp3) = exp3$Group.1  # 把分组列 Group.1（即 symbol）设为行名。
exp3 = exp3[,-c(1:4)]  # 删除前 4 列：Group.1（symbol）、Row.names（探针号）、probe_id、symbol，只保留样本表达列。
head(exp3)  # 打印前 6 行，人工核对行列结构是否正确。

# 整理分组信息
B <- cli[,c(1,2)]  # 提取样本信息的前两列；当前后续代码使用 title，并以样本行名对齐表达列，运行前应核对实际列名。
library(dplyr)  # 加载 dplyr，为下面的管道和 case_when() 准备。
B1 <- B %>%  # 用管道把左边的数据框传给右边函数的第一参数，%>% 来自 dplyr（magrittr），使代码自上而下连贯。
	mutate(group = case_when(  # mutate() 新建列；case_when() 是多条件分支赋值，类似 Excel 的多重 IF。
		grepl("NS", B$title) ~ "N",  # 若样本标题中包含字符串 NS，则 group 赋值为 "N"（正常对照）；grepl() 返回逻辑值，~ 左边是条件右边是结果。
		grepl("asthma", B$title) ~ "AS",  # 若标题包含 asthma，则赋值为 "AS"（哮喘组）；条件按顺序匹配，先命中者生效。
		TRUE ~ NA_character_  # 兜底规则：以上都不满足时赋 NA_character_（字符型缺失值），避免分组出错被悄悄归类。
	))  # 闭合 case_when() 与 mutate()。
group1 <- B1[order(B1$group),]  # 按 group 列排序（默认升序，会把 NA 排到最后），得到按分组排列的样本顺序表。
group1 <- group1[,-c(1:2),drop=FALSE]  # 删掉前两列（标题、样本号），只留 group 列；drop=FALSE 保证即使只剩一列也仍是数据框而不是向量。
exp4 <- exp3[,row.names(group1)]  # 按 group1 的行名（样本名）重排表达矩阵的列，使表达矩阵列顺序与分组表顺序完全一致——这是后续 PCA 分组着色与差异分析设计矩阵的前提。

# 保存表达矩阵
write.csv(exp4,"GSE179156_Matrix.csv")  # 把最终表达矩阵写成 CSV 保存到工作目录，默认带行名列（row.names=TRUE），作为本模块的交付结果。

# PCA主成分分析-----------------------------------------------------------------
gene <- exp4  # 把表达矩阵复制一份给 gene，避免直接修改原矩阵。
gene <- t(gene)  # 转置：PCA 要求「行 = 观测（样本）、列 = 变量（基因）」，转置后行为样本。

library(FactoMineR)  # 加载多元统计分析包，提供 PCA() 函数。
# 样本中基因表达值的 PCA 分析
gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)  # 执行 PCA：ncp=2 保留前两个主成分，scale.unit=TRUE 对各基因中心化并按标准差缩放，graph=FALSE 不自动绘图。
plot(gene.pca)  # 调用 PCA 结果的默认绘图方法，快速查看样本因子图；并不是单独绘制碎石图的命令。

# 提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])  # 从结果中取出每个样本在 Dim.1、Dim.2 上的坐标，转成数据框，这是画散点图的直接数据。
head(pca_sample)  # 查看前几行，确认坐标结构。

# 提取 PCA 前两轴的贡献度
pca_eig1 <- round(gene.pca$eig[1,2], 2)  # eig 矩阵第 2 列是「特征值占比百分比」，取第 1 主成分的百分比并四舍五入到 2 位小数。
pca_eig2 <- round(gene.pca$eig[2,2], 2)  # 同上，取第 2 主成分的百分比，用于标注坐标轴。

# 整理并合并样本分组信息
group <- group1[rownames(pca_sample), ]  # 按 PCA 样本行名提取分组，确保坐标和标签对应；group1 只剩一列时，默认降维会返回向量。
# group = c(rep("name1",5),rep("name2",2))  #自己构建分组信息
pca_sample <- cbind(pca_sample, group)  # 把分组列横向拼接到坐标数据框上，画图时才能按组着色。
pca_sample$samples <- rownames(pca_sample)  # 新增 samples 列保存样本名，供后面标签注释使用。
head(pca_sample)  # 查看合并后的数据：Dim.1、Dim.2、group、samples。

# ggplot2 绘制二维散点图
library(ggplot2)  # 加载绘图包。

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +  # 建立图层对象：横轴为第 1 主成分、纵轴为第 2 主成分，数据为 pca_sample。
	geom_point(aes(color = group), size = 3) +  # 画散点，颜色按 group 分组映射，点大小为 3。
	scale_color_manual(values = c("#ea5b57", "#fac03d", "#0c3c5f")) +  # 手动指定分组颜色（红、黄、深蓝三个十六进制色值），按因子水平顺序对应。
	theme(panel.grid = element_blank(),  # 主题设置开始：去掉网格线。
				panel.background = element_rect(color = 'black', fill = 'transparent'),  # 绘图区背景透明，并加黑色边框。
				legend.key = element_rect(fill = 'transparent')) +  # 图例色块背景透明，避免遮挡。
	labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  # 设置坐标轴标题：用 paste() 把「PCA1:」与贡献度百分比拼接；color='' 隐藏图例标题。
print(p)  # 在图形设备中渲染对象 p（交互窗口显示）。

# 将样本名称标记到图中
library(ggrepel)  # 加载标签防重叠包。

p1 <- p +  # 在原图 p 上叠加图层，结果存为 p1。
	geom_text_repel(aes(label = samples), size = 3, show.legend = FALSE,  # 添加样本名文本标签；geom_text_repel 会自动推开标签避免重叠；show.legend=FALSE 不为文字生成图例。
									box.padding = unit(0.5, 'lines'))  # 设置标签框周围的留白为 0.5 行高；box.padding 控制标签框留白，标签与数据点的留白另由 point.padding 控制。
print(p1)  # 显示带标签的版本。

# 添加 95% 置信椭圆
p + stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)  # 在 p 上叠加按组计算、level=0.95 的数据分布椭圆轮廓；该椭圆不是组均值的置信区间，也不是组间差异检验。
p_final <- p +  # 方式二（最终采用）：生成最终图对象 p_final。
	stat_ellipse(aes(fill = group), geom = "polygon", level = 0.95,  # 按组绘制 level=0.95 的分布椭圆，并用 polygon 进行填充；fill 按 group 映射。
							 alpha = 0.1, show.legend = FALSE) +  # 填充透明度 0.1（半透明，可透视下面的散点）；不生成图例。
	scale_fill_manual(values = c("#ea5b57", "#fac03d", "#0c3c5f"))  # 手动指定椭圆填充色，与散点颜色保持一致。

# 保存 PDF
pdf("PCA.pdf", width = 8, height = 6)  # 打开 PDF 图形设备，尺寸 8×6 英寸；之后画的图都写入该文件。
print(p_final)  # 把最终图输出到 PDF 设备。
dev.off()  # 关闭当前 PDF 设备，使文件正常写出；图形设备应与打开操作配对关闭。

# 保存高清 PNG
png("PCA.png", width = 2400, height = 1800, res = 300)  # 打开 2400×1800 像素、分辨率 300 dpi 的 PNG 设备，之后把图打印到这个设备。
print(p_final)  # 同一张图输出到 PNG 设备。
dev.off()  # 关闭 PNG 设备并保存。


# 构建分组矩阵
table(group1)  # 统计分组表中各组的样本数，用于核对下面硬编码的 57/29 是否正确。
group = c(rep("AS",57),rep("N",29))  # 手工构造 57 个 AS 和 29 个 N 的分组标签，必须逐列对应 exp4；样本筛选或顺序改变后要同步核对，不能只核对总数。
design = model.matrix(~0 + factor(group))  # 生成设计矩阵：~0+ 表示不含截距项（即不设默认参考组），每个组各占一列，1 表示该样本属于该组；limma 的 lmFit 必须用设计矩阵。
# model.matrix() 根据分组生成设计矩阵，用于线性模型拟合。
head(design)  # 查看前几行，重点检查列名。
colnames(design) = c("AS", "N")  # 按当前设计矩阵列顺序命名为 AS、N；先检查原列含义，避免仅改列名而把真实组别标反。
rownames(design) = colnames(exp4)  # 将设计矩阵行名设为表达矩阵列名，便于识别样本；改名字不会自动重新排列分组，顺序仍需预先核对。

# 进行差异表达分析
fit = lmFit(exp4, design)  # 对每个基因拟合线性模型：exp4 为表达矩阵（行=基因），design 为分组，得到每个基因的系数估计。
# design 为样本×模型系数的设计矩阵，编码分组及其他协变量。
cont.matrix = makeContrasts(contrasts="AS-N", levels=design)  # 构建对比向量 AS - N，即检验「AS 组减 N 组」的差异，正值代表在哮喘中高表达。
fit2 = contrasts.fit(fit, cont.matrix)  # 把线性模型按对比向量重新计算每个基因的对比系数与标准误，为检验做准备。

# 计算统计量和重要基因表
fit2 = eBayes(fit2)  # 经验贝叶斯 moderation：借用所有基因的信息收缩方差估计，提高小样本下检验的稳定性，是 limma 的核心步骤。
DEG_data1 = topTable(fit2, adjust="fdr", sort.by="B", number=nrow(exp4))  # 导出差异结果表：adjust="fdr" 用 BH 法校正 P 值（列 adj.P.Val）；sort.by="B" 按 B 统计量（log-odds）降序排；number=nrow(exp4) 返回全部基因而不是默认前 10 个。
#提取差异分析结果：sort.by="B"按照B统计量（也称为lods或log-odds）对基因进行排序

# 当前实际筛选：|logFC| >= 0.5 且原始 P.Value < 0.05；adj.P.Val 已计算，但未用于以下打标签。
DEG_data1$type = ifelse((DEG_data1$logFC >= 0.5 & DEG_data1$P.Value < 0.05),"up",  # 开始嵌套 ifelse()：logFC>=0.5 且原始 P.Value<0.05 标为 up。当前筛选使用原始 P 值，不是已计算的 adj.P.Val。
												ifelse((DEG_data1$logFC<= -0.5 & DEG_data1$P.Value < 0.05) ,"down","none"))  # 否则判断 logFC<=-0.5 且原始 P.Value<0.05，满足标为 down，其余标为 none。

DEG_data1$type1 = ifelse((DEG_data1$logFC >= 0.5 & DEG_data1$P.Value < 0.05),"up",  # 按相同阈值建立 type1，上调条件仍为 logFC>=0.5 且原始 P.Value<0.05；保留 up/down/none 原始标签供后续筛选。
												 ifelse((DEG_data1$logFC<= -0.5 & DEG_data1$P.Value < 0.05) ,"down","none"))  # 补全 type1 的下调和其余分支；后面 type 会改为带数量的因子标签，type1 则保留原始类别。

lab_title = c(paste0("up:",length(which(DEG_data1$type=="up"))),  # 统计上调基因数并拼成 "up:xx" 字符串。
							paste0("no change:",length(which(DEG_data1$type=="none"))),  # 统计不显著基因数，拼成 "no change:xx"。
							paste0("down:",length(which(DEG_data1$type=="down"))))  # 统计下调基因数，拼成 "down:xx"。这三个字符串将作为火山图图例标签。

DEG_data1$type = factor(DEG_data1$type,  # 把 type 列转成因子并自定义水平顺序。
												levels = c('up','none','down'),  # 固定水平顺序为 up、none、down，使图例与颜色按该顺序出现。
												labels = lab_title)  # 把水平名替换为含数量的标签（up:xx / no change:xx / down:xx）。

write.csv(DEG_data1,"DEG_GSE179156.csv")  # 保存完整差异基因表（含 logFC、AveExpr、t、P.Value、adj.P.Val、B、type、type1）。


# 差异基因火山图----------------------------------------------------------------
p_DEG = ggplot(DEG_data1, aes(logFC, -log10(P.Value),color = type))  # 预先定义一个图层对象（美学映射：x=logFC，y=−log10P，颜色=分组标签）；本行定义的 p_DEG 在下文未被使用。

p<- ggplot(DEG_data1,aes(x=logFC,y=-log10(P.Value)))+  # 正式绘图对象 p：横轴 logFC，纵轴为 P 值取 −log10（值越大越显著）。
	scale_x_continuous(limits=c(-3,4), breaks=seq(-3,4,1))+  # 限定横轴范围 −3 到 4，刻度间隔 1（seq(-3,4,1) 生成等差序列）；超出范围的点会被去掉。
	geom_point(aes(color=type),size =1.5)+  # 画点，颜色按 type（up/none/down）映射，点大小 1.5。
	scale_color_manual(values=c("red","gray","blue"),labels = c("Up", "none", "Down"))+  # 手动指定颜色与图例文字：按因子顺序 up=红、none=灰、down=蓝。
	geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8,col="#4C5B61")+  # 画水平参考线于 −log10(0.05)≈1.30（显著性阈值）；lty=4 虚线、lwd 线宽、alpha 透明度。
	geom_vline(xintercept = c(0.5,-0.5),lty=4,lwd=0.6,alpha=0.8,col="#4C5B61")+  # 在 logFC = ±0.5 处画两条竖线，即差异倍数阈值，与打标签的标准一致。
	annotate("text",x=2,y=11,label=lab_title[1],size=5,col="black") +  # 在图内 (2, 11) 位置添加文本，内容为 lab_title[1]（上调基因数）。
	annotate("text",x=-2,y=11,label=lab_title[3],size=5,col="black") +  # 在左侧 (−2, 11) 添加文本，内容为 lab_title[3]（下调基因数）。
	theme_bw()+  # 使用黑白主题（白底+边框）。
	theme(panel.grid=element_blank(),axis.text = element_text(size = 12),	axis.title = element_text(size = 12),legend.text = element_text(size = 12))  # 去掉网格线，并把坐标轴刻度文字、轴标题、图例文字统一设为 12 号字。
png("vlocano.png",width = 7,height = 6,units = "in",res = 300)  # 打开 7×6 英寸、300 dpi 的 PNG 设备；当前下一行立即打开 PDF，PNG 在本段没有得到绘制内容。
pdf("vlocano.pdf",width = 7,height = 6)  # 打开 PDF 并成为当前设备，前面的 PNG 转到后台；正确导出两种格式时，每个设备都应分别打开、print(p)、dev.off()。
p  # 交互运行时显示 p，当前设备为 PDF；整份脚本用 source() 运行时应显式 print(p)，单独写对象名通常不会自动绘制。
dev.off()  # 关闭当前 PDF 设备；本段先打开的 PNG 设备仍未配对关闭，也尚未绘制。

# 聚类分析热图------------------------------------------------------------------
#BiocManager::install("ggvenn")
library(ComplexHeatmap)  # 加载 ComplexHeatmap 包，功能强大的热图绘制工具，支持注释条、聚类等。
library(pheatmap)  # 加载 pheatmap（本脚本实际使用的是 ComplexHeatmap，此行备用）。

# 提取全部上下调差异基因的表达矩阵。
down <- DEG_data1[DEG_data1$type1 == 'down',]  # 用 type1 筛出全部下调基因行（用备份列而不用被改写标签的 type）。
up <- DEG_data1[DEG_data1$type1 == 'up',]  # 筛出全部上调基因行。
up_down <- rbind(up,down)  # 上下调基因表纵向合并，得到「全部差异基因」表；这个 up_down 后面 GO/KEGG 富集也用它。
exp_updown <- exp4[rownames(up_down),]  # 按差异基因名从表达矩阵取子矩阵（全部差异基因的表达）。

# 按 logFC 大小选择上调、下调各 20 个基因的表达矩阵。
df_2 <- DEG_data1[order(DEG_data1$logFC),]  # 按 logFC 升序排列：最前面是最负（下调最强），最后面是最正（上调最强）。
down_20gene <- df_2[grep("down",df_2$type1)[1:20],]  # 在已按 logFC 升序排列的下调基因中取前 20 个；这里按变化幅度而不是 P 值排序，若不足 20 个，固定下标可能引入 NA 行。
up_20gene <- tail(df_2[grep("up",df_2$type1),], n =20)  # 从上调基因中取 logFC 最大的末尾 20 个；若上调基因不足 20 个，tail() 只返回现有行。
updown_20gene <- rbind(up_20gene,down_20gene)  # 合并上调和下调所选基因；当两类均足够时合计 40 个，而不是一共 20 个。
exp_updown_20gene <- exp4[rownames(updown_20gene),]  # 提取这 40 个基因在全部样本中的表达子矩阵。

C <- exp_updown_20gene  # 赋值给 C，作为热图输入矩阵。
C <- as.matrix(C)  # 转成矩阵类型，scale() 与 Heatmap 需要矩阵。
samples <- group  # 取分组向量，后面作为顶部注释条。
for (i in 1:nrow(C)) C[i, ] <- scale(as.numeric(C[i, ]))  # 逐行做 z-score 标准化：减去该基因均值除以标准差，使不同表达量的基因能在同一色阶上比较（热图标配步骤）。
C[is.na(C)]=0  # 把矩阵中的 NA/NaN 置为 0，其中可能包括零方差基因标准化产生的无效值；原有缺失也会一起被替换。

D <- Heatmap(C,  # 调用 Heatmap() 开始构建热图对象 D，第一个参数是矩阵。
						 name =  "expression",  # 图例名称设为 “expression”。
						 col = colorRampPalette(c("navy","white","firebrick3"))(100),  # 用蓝—白—红生成 100 个渐变色：低表达蓝、中间白、高表达红（发散配色，适合 z-score）。
						 show_row_names = T,  # T 表示 TRUE，因此显示基因行名；原先“不展示行名”的注释写反了。
						 show_column_names = F,  # 不显示列名（样本名太多，避免拥挤）。
						 column_names_gp = gpar(fontsize = 8),  # 列名（样本）字体 8 号；gpar() 来自 grid。
						 column_names_rot = 90,  # 列名旋转 90 度竖排。
						 use_raster = F,  # 不主动将热图主体栅格化；PDF 中通常保留矢量元素，PNG 最终仍是位图。TRUE 可减小大型热图的矢量输出体积。
             #cluster_columns = F,
						 row_names_gp = gpar(fontsize = 8),  # 行名（基因名）字体 8 号。
						 top_annotation = HeatmapAnnotation(Group = samples,  # 顶部添加注释条，变量名 Group，取值为分组向量。
						 																	 simple_anno_size = unit(2, 'mm'),  # 注释条高度为 2 毫米。
						 																	 col = list(Group = c('N' = '#00DAE0',  # 为注释条指定颜色映射：N 组为青色。
						 																	 										 'AS' = '#FF9289')),  # AS 组为橙红色。
						 																	 show_annotation_name = FALSE))  # 不显示注释条左侧的名字，收尾。
D  # 交互运行时显示热图 D；在 source() 或非交互脚本中，通常应显式使用 ComplexHeatmap::draw(D)。
pdf("heatmap_20.pdf",width = 7,height = 6)  # 打开热图 PDF 设备；当前下一行紧接着打开 PNG，PDF 尚未接收到热图。
png("heatmap_20.png",width = 7,height = 6,units = "in",res = 300)  # 打开 PNG 并切为当前设备；要导出两种格式，应该为每个设备分别打开、draw(D)、dev.off()。
D  # 交互运行时在当前 PNG 设备显示 D；整份脚本运行时应显式 draw(D)，避免只写对象名而未绘图。
dev.off()  # 关闭当前 PNG 设备；前面打开的热图 PDF 仍未配对关闭，也未在该设备绘制。


# Venn 维恩图-------------------------------------------------------------------
# BiocManager::install(c('ggvenn','VennDiagram'))
library(ggvenn)  # 加载 ggvenn，绘制出版风格韦恩图（本脚本绘图用它）。
library(VennDiagram)  # 加载 VennDiagram，本脚本用它的 get.venn.partitions() 提取交集成员。

# 读取数据文件
venn_dat  <- read.csv("~/新转录组基础班20260724/文章复现2026724/8. venn图/Gene.csv",header = T)  # 读取课程目录中的 Gene.csv；~ 是当前用户主目录，header=T 表示首行为集合名称。换环境需核对实际路径及各列基因来源。
venn_list <- as.list(venn_dat)  # 把数据框转成列表，每个元素是一个集合的基因向量，这是 ggvenn() 需要的输入格式。
venn_list <- purrr::map(venn_list, na.omit)  # 用 purrr 的 map() 对列表每个元素执行 na.omit()，删除缺失值。
venn_list <- purrr::map(venn_list, function(x){x[x!=""]})  # 对列表中每个基因向量去掉空字符串；这一步没有清除空格或统一基因命名，也不是执行去重。

# 绘图
ggvenn(  # 开始调用 ggvenn 绘图。
	data = venn_list,  # 传入整理好的基因集合列表。
  #columns = NULL,           # 对选中的列名绘图，最多选择4个，NULL为默认全选  # 备选列选择参数，当前已注释；NULL 通常使用全部列，可用列名指定需要展示的集合。
	show_elements = F,  # FALSE 时区域里显示元素个数而非基因名。
	label_sep = "\n",  # 基因名之间的分隔符为换行（仅 show_elements=TRUE 时生效）。
	show_percentage = F,  # 不显示各区域占比百分比。
	digits = 1,  # 设置百分比标签的小数位数为 1；当前 show_percentage=F，因此不显示百分比。
	fill_color = c("#E41A1C", "#1E90FF", "#FF8C00", "#80FF00"),  # 提供可用的集合填充色；当前真正使用的颜色数量取决于参与绘图的集合数。
	fill_alpha = 0.5,  # 半透明填充，交集处颜色会叠加。
	stroke_color = "white",  # 边框白色。
	stroke_alpha = 0.5,  # 边框半透明。
	stroke_size = 0.5,  # 边框线宽 0.5。
	stroke_linetype = "solid",  # 实线边框；注释列出可选线型（twodash、dotted 等）。
	set_name_color = "black",  # 集合名称黑色。
	set_name_size = 6,  # 集合名称字号 6。
	text_color = "black",  # 数字标签黑色。
	text_size = 4  # 数字字号 4。
)  # 结束 ggvenn() 调用；交互运行时显示图。当前没有显式保存图对象，整份脚本运行后再用 ggsave() 的行为应核对。

# 保存图片
ggsave("venn.png", width = 6, height = 6, dpi = 300,bg = 'white')  # 把上一个 ggplot 对象存为 6×6 英寸、300 dpi 的 PNG，白色背景。
ggsave("venn.pdf", width = 6, height = 6)  # 同一张图存为矢量 PDF。

#提取组间交集元素
inter <- get.venn.partitions(venn_list)  # 计算韦恩图每个分区的成员信息，返回数据框：包含各分区的集合组合、基因列表、数量等。
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')  # 循环遍历每一行，把 ..values.. 列中的基因向量用逗号拼成一个字符串，写入 values 列（方便 Excel 查看）。
#目的是遍历数据框inter的每一行，将每行中列名为..values..的那一列数据拼接成一个以逗号分隔的字符串，
#然后赋值给同一行中列名为values的单元格
write.csv(inter[-c(3, 4)], 'venn_inter.csv', row.names = FALSE)  # 删除 inter 的第 3、4 列后导出，不写行名；列位置的具体含义取决于集合数和返回结构，应先检查 names(inter)。

#GO富集分析------------------------------------------------------------
# 加载富集分析所需要的R包
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)  # 加载富集分析核心包，提供 enrichGO()、enrichKEGG()、bitr()。
library(org.Hs.eg.db)  # 加载人类注释数据库，支持 SYMBOL、ENTREZID 等标识之间的映射。
library(ggplot2)  # 加载绘图包（富集结果图基于 ggplot2）。

Genes <- bitr(rownames(up_down),  # 把全部上下调差异基因的 SYMBOL 转为 ENTREZID，不是只转换热图或 Venn 交集基因；需留意未映射和一对多映射。
							fromType = "SYMBOL",  # 声明输入 ID 类型是基因 symbol。
							toType = c("ENTREZID"),  # 转换成 ENTREZID（富集分析数据库用的编号）。
							OrgDb = org.Hs.eg.db)  # 指定注释库为人类，保证转换依据正确。

# GO富集分析与结果存储
GO <- enrichGO(gene = Genes$ENTREZID,  # 用转换后的差异基因做 GO 过度富集分析；当前未传入 universe，背景使用数据库默认的可注释基因范围。
							 OrgDb = org.Hs.eg.db,  # 指定注释数据库。
							 keyType = "ENTREZID",  # 声明输入基因 ID 类型。
							 ont = "all",  # 同时分析 BP、CC、MF 三个本体；也可只填其中一个。
							 pAdjustMethod = "BH",  # 用 Benjamini-Hochberg 法校正多重检验。
							 pvalueCutoff = 1,  # 将 P 值相关过滤阈值设为 1，减少阈值过滤；保留在结果中的条目仍需按校正 P 值判断是否显著。
							 qvalueCutoff = 1,  # 将 q 值过滤阈值设为 1；这不等于每个 GO 条目都会返回，结果仍受映射、重叠和集合大小等条件限制。
							 minGSSize = 5,  # 要求背景范围中的 GO 基因集合至少包含 5 个基因；限制的是条目集合大小，不是命中的差异基因至少有 5 个。
							 maxGSSize = 5000,  # 要求背景范围中的 GO 基因集合最多包含 5000 个基因，与 minGSSize 一起限定参与检验的集合大小。
							 readable = TRUE)  # 结果中自动把 ENTREZID 转成基因 symbol，便于阅读。
write.csv(data.frame(ID=row.names(GO@result),GO@result),  # 读取富集结果的 @result 表并另加一列行名 ID 后导出；@result 本身通常已有 ID 列，因此可能出现重复 ID 列。
					file = "GO_result_GSE179156.csv")  # 保存为 CSV（含 ID、Description、GeneRatio、pvalue、p.adjust、qvalue、基因列表等列）。

# 结果可视化--柱状图
pdf(file="GO_bar_GSE179156.pdf",width = 10,height = 10)  # 打开 PDF 设备，10×10 英寸。
# png("GO柱状图_LJ_2.png",width = 10,height = 10,units = "in",res = 300)
barplot(GO, drop = TRUE,  # 绘制 GO 富集柱状图；条目通常在纵轴，柱长使用默认的 Count，具体以图轴标签和所用版本为准，不能直接说成 GeneRatio。
				showCategory =10,  # 每个拆分组最多展示 10 个条目；这是展示数量设置，不是显著性阈值，返回不足时按实际数量显示。
				label_format=100,  # 条目名称文本最长 100 字符后换行，防止标签过长溢出。
				split="ONTOLOGY") +  # 按 ONTOLOGY 列（BP/CC/MF）拆分成多个面板。
	facet_grid(ONTOLOGY~., scale='free') +  # 按 ONTOLOGY 纵向分面，scale 部分匹配 scales；free 允许坐标尺度变化，但 facet_grid 同列仍共享横轴尺度。
	theme(axis.text.x = element_text(size = 10),  # x 轴刻度字号 10。
				axis.text.y = element_text(size = 10),  # y 轴刻度字号 10。
				strip.text = element_text(size = 8))  # 分面标签（BP/CC/MF）字号 8。
#panel.spacing = unit(1, "lines")) #调整面板之间的间距
#axis.title.x = element_text(size = 14),  # 调整x轴标题字体大小
#axis.title.y = element_text(size = 14),  # 调整y轴标题字体大小
#legend.title = element_text(size = 12),  # 调整图例标题字体大小
#legend.text = element_text(size = 10))   # 调整图例文本字体大小
dev.off()  # 关闭 PDF 设备，保存图片。
head(GO@result)  # 查看富集结果前 6 行，检查条目与 P 值。

# 结果可视化--点图
pdf(file="GO_dot_GSE179156.pdf",width = 10,height = 10)  # 打开 PDF 设备。
# 方式一：分面
dotplot(GO,showCategory = 10,label_format=100,split="ONTOLOGY") +  # 绘制 GO 点图：默认横轴 GeneRatio、点大小 Count、颜色 p.adjust；按本体拆分，每类最多显示 10 项。展示条目不等于都显著。
	facet_grid(ONTOLOGY~., scale='free') +  # 纵向分面 BP/CC/MF。
	scale_size(range = c(2, 4)) +  # 点大小映射范围缩放到 2–4，使视觉差异更温和。
	theme(  # 开始主题设置。
    #text = element_text(size = 14),  # 修改整体字体大小
    #axis.title = element_text(size = 16),  # 修改轴标题字体大小
		axis.text.y = element_text(size = 10),  # y 轴（条目名）字号 10。
		axis.text.x = element_text(size = 10)  # x 轴字号 10。
    #legend.title = element_text(size = 14),  # 修改图例标题字体大小
    #legend.text = element_text(size = 12)   # 修改图例文本字体大小
	)  # 闭合 theme()。
dev.off()  # 关闭设备，保存点图。

pdf(file="GO_GSE179156_1.pdf",width = 10,height = 10)  # 打开 PDF 设备。
# 方式二：不分面
dotplot(  # 开始绘制不分面的点图。
	GO,  # 富集结果对象。
	x = "GeneRatio",  # 横轴为 GeneRatio，即命中该条目的输入基因数/该分析中可注释的输入基因总数；分母不是该条目的所有基因数。
	color = "p.adjust",  # 按校正后 P 值 p.adjust 映射颜色；数值与颜色的对应关系以实际图例为准。
	showCategory = 15,  # 最多展示 15 个条目；该参数控制显示数量，不额外限定统计显著性。
	size = NULL,  # 点大小用默认映射（基因数）。
	split = NULL,  # 不按结果列拆分条目；本段也未添加 facet_grid，因此显示为单一面板。
	font.size = 12,  # 字号 12。
	title = "",  # 不设图标题。
	orderBy = "x",  # 按 x 轴值（GeneRatio）排序展示。
	label_format = 30)  # 条目名超过 30 字符换行。
dev.off()  # 关闭设备，保存 PDF。

# KEGG富集分析与结果存储 -------------------------------------------------------
KEGG <- enrichKEGG(gene = Genes$ENTREZID,  # 执行 KEGG 通路富集，输入与 GO 相同的差异基因 ENTREZID。
									 organism = "hsa",  # 物种代码：hsa 为人；小鼠为 mmu（注释给出提示）。
									 keyType = "kegg",  # 按 KEGG ID 类型读取；在人类 hsa 中通常使用对应的 ENTREZID 数字编号，其他物种不能一概照搬。
									 pAdjustMethod = "BH",  # BH 法校正多重检验。
									 pvalueCutoff = 1,  # 把 KEGG 的 P 值相关过滤阈值设为 1，保留较完整结果；仍需检查 p.adjust，不能把所有返回通路称为显著富集。
									 qvalueCutoff = 1)  # 把 q 值过滤阈值设为 1 并结束调用；未提供 universe 时使用数据库默认背景，在线注释来源和版本也会影响结果。

barplot(KEGG,  # 绘制 KEGG 柱状图。
				x = "GeneRatio",  # 横轴选择 GeneRatio：该通路命中的输入基因数/参与该分析的可注释输入基因总数。
				label_format =  100,  # 通路名最长 100 字符后换行。
				color = "p.adjust",  # 柱色按校正 P 值渐变。
				showCategory = 15,  # 展示前 15 条通路。
				title = "KEGG_enrichment")  # 设置图形标题为 KEGG_enrichment。
dotplot(KEGG,  # 绘制 KEGG 点图；本段没有为 KEGG 图单独打开文件设备或调用 ggsave()，因此并未明确导出 KEGG 图片文件。
				label_format =  100  # 通路名换行阈值 100 字符。
)  # 结束点图调用（此图输出到当前设备）。

# 结果存储
KEGG_results <- clusterProfiler::setReadable(KEGG,  # 使用当前代码的 clusterProfiler::setReadable() 把 KEGG 结果中的基因编号转换为可读基因符号；保留当前命名空间写法。
																	OrgDb = org.Hs.eg.db,  # 指定人类注释库。
																	keyType = "ENTREZID")  # 声明要转换的 ID 类型。
write.csv(data.frame(ID=rownames(KEGG_results@result),KEGG_results@result),  # 把可读化的 KEGG @result 导出并另加行名 ID；结果表通常已有 ID 列，需留意重复 ID 列。
					file = "KEGG_enrich_result_GSE179156.csv")  # 保存 KEGG 富集结果表，脚本到此结束。

