# 课程讲解版：蛋白组代码讲解7.3
# 模块开头说明目的；行内注释逐个说明函数作用，并保留本行用途和关键设置。
# 常用符号：<- 或 = 保存结果；[] 取子集；$ 取字段；%>% 将数据传给下一步；ggplot 中 + 叠加图层。

# 加载包

# ---- 01. 加载基础工具 ----
# 目的：准备数据整理、统计分析和绘图所需软件包。
library(GEOquery)  # library() 加载 GEOquery，用于读取 GEO 数据。
library(stringr)  # library() 加载 stringr，用于处理字符串。
library(ggplot2)  # library() 加载 ggplot2，用于按图层绘图。
library(dplyr)  # library() 加载 dplyr，用于整理、筛选和分组数据。
library(ggpubr)  # library() 加载 ggpubr，用于绘制统计图并添加比较结果。
library(grid)  # library() 加载 grid，用于设置图形尺寸和字体属性。
library(ggrepel)  # library() 加载 ggrepel，用于添加自动避让的文字标签。
library(limma)  # library() 加载 limma，用于进行线性模型差异分析。

# PCA主成分分析-----------------------------------------------------------------
# 读取蛋白表达值矩阵

# ---- 02. PCA 输入整理 ----
# 目的：汇总同名基因并转为样本×蛋白，为样本层面 PCA 准备数据。
setwd("~/蛋白组课程-6.17")  # setwd() 将工作目录设为蛋白组课程文件夹，后续使用相对路径读写。
pro <- read.csv('1.PCA分析/可信蛋白.csv',header = T)  # read.csv() 读取可信蛋白定量表，header=TRUE 表示首行为列名。
pro1 = aggregate(x=pro,by=list(pro$Gene.Name),FUN = max)  # aggregate() 按 Gene.Name 分组，对同名记录各列取最大值；max() 取最大值；list() 将对象组成列表。
rownames(pro1) = pro1$Group.1  # rownames() 读取或设置行名；把分组后的基因符号设为行名。
pro1 = pro1[,-c(1:2)]  # c() 将多个值合成向量；删除前两列分组和注释信息，保留样本表达列。
pro1 <- t(pro1)  # t() 转置行列；转置为样本×蛋白，供 PCA 分析样本差异。


# ---- 03. 计算 PCA ----
# 目的：用前两个主成分概括整体表达差异，查看样本结构。
library(FactoMineR)  # library() 加载 FactoMineR，用于计算 PCA 等多元分析。
# 样本中蛋白表达值的 PCA 分析
pro.pca <- PCA(pro1, ncp = 2, scale.unit = TRUE, graph = FALSE)  # PCA() 计算前两个主成分；按蛋白标准化，关闭自动绘图。
plot(pro.pca)  # plot() 显示 PCA 样本坐标的简图。

# 提取样本在 PCA 前两轴中的坐标

# ---- 04. 提取坐标和解释率 ----
# 目的：整理样本位置及两轴方差解释率，供绘图和轴标题使用。
pca_sample <- data.frame(pro.pca$ind$coord[ ,1:2])  # data.frame() 创建数据框；提取各样本在 PC1、PC2 上的坐标，转为数据框。
head(pca_sample)  # head() 查看数据前六行，核对结构和内容。

# 提取 PCA 前两轴的贡献度
pca_eig1 <- round(pro.pca$eig[1,2], 2)  # round() 提取 PC1 的方差解释率，保留两位小数。
pca_eig2 <- round(pro.pca$eig[2,2], 2)  # round() 提取 PC2 的方差解释率，保留两位小数。

# 整理并合并样本分组信息
# group <- group[rownames(pca_sample), ]  # 备用代码，未执行；rownames() 读取或设置行名。

# ---- 05. 对齐分组标签 ----
# 目的：将样本坐标与实际组别对应，准备按组着色。
group = c(rep("H_a",4),rep("L_a",4))  # c() 将多个值合成向量；rep() 重复指定内容；构造前四个 H_a、后四个 L_a 标签，须与样本顺序一致。
pca_sample <- cbind(pca_sample, group)  # cbind() 将分组信息按列拼接到 PCA 坐标表。
pca_sample$samples <- rownames(pca_sample)  # rownames() 读取或设置行名；新增 samples 列保存样本名，供绘图标注。
head(pca_sample)  # head() 查看数据前六行，核对结构和内容。

# ggplot2 绘制二维散点图

# ---- 06. PCA 基础散点图 ----
# 目的：展示两组样本在前两个主成分上的位置。
library(ggplot2)  # library() 加载 ggplot2，用于按图层绘图。

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +  # ggplot() 建立 PCA 散点图，横轴为 PC1，纵轴为 PC2；aes() 将数据变量映射到坐标或外观。
  geom_point(aes(color = group), size = 3) +  # geom_point() 每个点代表一个样本，按 group 着色，点大小设为 3；aes() 将数据变量映射到坐标或外观。
  scale_color_manual(values = c('red', 'green','blue')) +  # scale_color_manual() 手动指定分组颜色；当前两组只需使用其中两色；c() 将多个值合成向量。
  theme(panel.grid = element_blank(),  # theme() 去掉绘图区网格线，开始设置主题；element_blank() 隐藏相应主题元素。
        panel.background = element_rect(color = 'black', fill = 'transparent'),  # element_rect() 设置透明背景和黑色绘图区边框。
        legend.key = element_rect(fill = 'transparent')) +  # element_rect() 将图例色块背景设为透明。
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  # labs() 轴标题加入 PCA 解释率，并隐藏颜色图例标题；paste() 按分隔符拼接文字。
print(p)  # print() 显示图形，并输出到当前图形设备。

# 将样本名称标记到图中

# ---- 07. PCA 样本标签 ----
# 目的：为各点添加样本名称，便于定位聚集或偏离的样本。
library(ggrepel)  # library() 加载 ggrepel，用于添加自动避让的文字标签。

p1 <- p +  # 在原 PCA 图 p 上继续添加标签图层。
  geom_text_repel(aes(label = samples), size = 5, show.legend = FALSE,  # geom_text_repel() 添加自动避让的样本名标签，字号为 5，隐藏文字图例；aes() 将数据变量映射到坐标或外观。
                  box.padding = unit(0.5, 'lines'))  # unit() 给长度附加单位；设置标签框四周留白为 0.5 行高。
print(p1)  # print() 显示图形，并输出到当前图形设备。

# 添加 level=0.95 的数据分布椭圆；不是组均值置信区间，也不是组间检验。

# ---- 08. PCA 分布椭圆 ----
# 目的：用组内分布椭圆辅助观察样本散布范围。
p + stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)  # stat_ellipse() 按组绘制数据分布椭圆；aes() 将数据变量映射到坐标或外观；按组添加 level=0.95 的分布椭圆，不作组间显著性判断。
p1 + stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +  # stat_ellipse() 按组绘制数据分布椭圆；aes() 将数据变量映射到坐标或外观；添加填充分布椭圆，透明度为 0.1，隐藏该层图例。
  scale_fill_manual(values = c('red', 'green','blue'))  # scale_fill_manual() 手动指定椭圆填充色，与样本分组对应；c() 将多个值合成向量。
#stat_ellipse可以根据数据的均值和协方差矩阵计算出椭圆的参数,并将其绘制在散点图或密度图上.通过调整参数,
#可以控制椭圆的大小和形状，以及绘制线条的颜色和风格

# 样品间相关性分析####
# 相关性分析
#BiocManager::install("GGally")  # 备用代码，未执行；BiocManager::install() 安装或更新软件包。

# ---- 09. 相关分析输入 ----
# 目的：读取独立蛋白表并去重，准备计算样本之间的相关。
library(corrplot)  # library() 加载 corrplot，用于绘制相关矩阵。
library(ggpubr)  # library() 加载 ggpubr，用于绘制统计图并添加比较结果。
library(ggplot2)  # library() 加载 ggplot2，用于按图层绘图。
library(circlize)  # library() 加载 circlize，用于进行环形绘图及状态管理。
# 读取txt数据
tf <- read.csv("2.相关性分析/相关性分析.csv",header = T)  # read.csv() 读取相关性分析的蛋白表达表，首行为列名。
tf <- tf[!duplicated(tf$Gene.Name),]  # duplicated() 标记重复记录；按 Gene.Name 去重，同名记录只保留第一条。
rownames(tf) = tf$Gene.Name  # rownames() 读取或设置行名；将 Gene.Name 设为行名。
tf=tf[,-1]  # 删除标识列，保留用于相关计算的样本表达列。


# ---- 10. 计算与预览样本相关 ----
# 目的：用 Pearson 系数概括样本表达模式的相近程度。
data <- cor(tf, method="pearson",use = "pairwise.complete.obs")  # cor() 计算样本间 Pearson 相关，每对样本使用共同非缺失值。
corrplot(data)  # corrplot() 按默认样式显示样本相关矩阵。

corrplot(data, method = "circle",  # corrplot() 用圆形展示样本相关系数。
         tl.col = "black", tl.cex = 1.2, tl.srt = 45)  # 标签设为黑色，缩放 1.2 倍，旋转 45 度。
corrplot(data, method = "ellipse",  # corrplot() 用椭圆形状展示相关强度与方向。
         type = "upper",  # 只绘制相关矩阵的上三角。
         tl.col = "black", tl.cex = 1.2, tl.srt = 45  # 标签设为黑色，缩放 1.2 倍，旋转 45 度。
)  # 结束本次函数调用。
# method和type还有其他选项，根据需要修改
#method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")  # 备用代码，未执行；c() 将多个值合成向量。
#type = c("full", "lower", "upper")  # 备用代码，未执行；c() 将多个值合成向量。


# ---- 11. 相关检验和图形输出 ----
# 目的：为相关矩阵设置色阶与检验标记，并组合上下三角显示。
col1 <- colorRampPalette(c("blue", "white", "red"))(100)  # colorRampPalette() 创建渐变配色函数；c() 将多个值合成向量；生成蓝—白—红的 100 色渐变。
#col2 <- colorRampPalette(c("white", "steelblue"))(100) # 单色渐变（适合显著性矩阵）  # 备用代码，未执行；colorRampPalette() 创建渐变配色函数；c() 将多个值合成向量。
#col3 <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(100) # 离散色阶（适合分类数据）  # 备用代码，未执行；colorRampPalette() 创建渐变配色函数；RColorBrewer::brewer.pal() 从命名调色板提取颜色。

testRes = cor.mtest(tf, method="spearman")  # cor.mtest() 计算 Spearman 检验 P 值；与前面的 Pearson 系数方法不同。

png("corrplot.png",width = 7,height = 6,units = "in",res = 300)  # png() 打开 7×6 英寸、300 dpi 的 PNG 图形设备。
pdf("corrplot.pdf",width = 7,height = 6)  # pdf() 打开 7×6 英寸 PDF 设备；此前 PNG 设备仍未关闭。
par(mar = c(5,2,2,2) + 0.3,  # par() 设置绘图区下、左、上、右边距；c() 将多个值合成向量。
    oma = c(2,1,1,1),  # c() 将多个值合成向量；设置图形外边距为 2、1、1、1 行。
    xpd = NA)  # 放宽绘图区裁剪，使文字可画到设备边界。
corrplot(data, method = "circle", col = col1,  # corrplot() 用自定义色阶绘制圆形相关图。
         tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt",  # 标签放在左侧和顶部，缩放 0.8 倍并旋转 45 度。
         p.mat = testRes$p, diag = T, type = 'upper',  # 传入检验 P 值，显示上三角及对角线。
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,  # c() 将多个值合成向量；按 0.001、0.01、0.05 分级标记显著性，符号缩放 1.2 倍。
         insig = 'label_sig', pch.col = 'grey20')  # 用深灰色符号标注相关检验的显著性等级。
corrplot(data, method = "number", type = "lower",col = col1,  # corrplot() 在下三角用数字展示同一相关矩阵。
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",  # 隐藏数字层的行列标签，避免重复显示。
         add = T)  # 将下三角数字叠加到已有相关图上。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。
circos.clear()  # circos.clear() 清除 circlize 绘图状态，便于后续重新绘图。

# 样品层次聚类分析####

# ---- 12. 样本聚类输入 ----
# 目的：整理样本×蛋白表达，并设置本段使用的工具和选项。
library(WGCNA)  # library() 加载 WGCNA，用于进行样本检查、聚类和共表达分析。
options(stringsAsFactors = FALSE)  # options() 关闭字符列自动转因子；此项并非开启多线程。
protein = read.csv("3.样品层次聚类分析/可信蛋白.csv",header = T)  # read.csv() 读取样本层次聚类所用的蛋白定量表。
protein <- protein[!duplicated(protein$protein.Name),]  # duplicated() 标记重复记录；按 protein.Name 去重，只保留首次出现的记录。
rownames(protein) = protein$protein.Name  # rownames() 读取或设置行名；将 protein.Name 设为行名。
protein1 = as.data.frame(t(protein[,-1]))  # as.data.frame() 转为数据框；t() 转置行列；去掉标识列并转置，得到样本×蛋白的数据框。
dim(protein1)  # dim() 查看数据行数和列数，核对输入规模。

# 数据处理
# 检查缺失值过多的蛋白和样本

# ---- 13. 过滤无效输入 ----
# 目的：删除缺失或变异检查未通过的蛋白与样本。
gsg = goodSamplesGenes(protein1, verbose = 3)  # goodSamplesGenes() 检查缺失值过多及变异不足的蛋白和样本，输出检查信息。
gsg$allOK  # 查看是否全部通过检查；TRUE 表示无需按本次结果过滤。

# 去除有问题的蛋白和样本
if (!gsg$allOK){  # 若存在未通过检查的样本或蛋白，则进入过滤步骤。
  if (sum(!gsg$goodGenes)>0)  # sum() 求和，逻辑值 TRUE 计为 1；判断是否存在不合格蛋白。
    printFlush(paste("Removing genes:", paste(names(protein1)[!gsg$goodGenes], collapse = ", ")))  # printFlush() 打印将被删除的蛋白名称，用逗号连接；paste() 按分隔符拼接文字；names() 读取或设置对象名称。
  if (sum(!gsg$goodSamples)>0)  # sum() 求和，逻辑值 TRUE 计为 1；判断是否存在不合格样本。
    printFlush(paste("Removing samples:", paste(rownames(protein1)[!gsg$goodSamples], collapse = ", ")))  # printFlush() 打印将被删除的样本名称；paste() 按分隔符拼接文字；rownames() 读取或设置行名。
  protein1 = protein1[gsg$goodSamples, gsg$goodGenes]  # 同时保留合格的样本行和蛋白列。
}  # 结束输入过滤的条件分支。


# ---- 14. 标准化与样本树 ----
# 目的：按蛋白标准化后计算样本距离，建立层次聚类树。
protein1.scaled <- scale(protein1)  # scale() 默认按列减均值并除以标准差；对每个蛋白跨样本做 z-score 标准化。

h.clust<-hclust(dist(protein1.scaled,method="euclidean"))  # hclust() 按欧氏距离建样本树，默认使用 complete 连接法；dist() 计算各行对象之间的距离。
sampleTree = hclust(dist(protein1.scaled), method = "average")  # hclust() 按欧氏距离和 average 连接法另建样本树；dist() 计算各行对象之间的距离。
# method="average\euclidean")

# ---- 15. 展示样本聚类 ----
# 目的：查看样本间的层次关系及可能偏离主簇的样本。
png("hclust.png",width = 7,height = 6,units = "in",res = 300)  # png() 打开 7×6 英寸、300 dpi 的 PNG 图形设备。
pdf("hclust.pdf",width = 7,height = 6)  # pdf() 打开 7×6 英寸的 PDF 图形设备。
plot(sampleTree,hang=-1,main="Cluster Analysis Of Samples",  # plot() 绘制样本聚类树，hang=-1 使叶端对齐。
     sub="",  # 隐藏聚类图副标题。
     xlab="",  # 隐藏横轴标题。
     ylim = c(0, max(h.clust$height)+1),  # c() 将多个值合成向量；max() 取最大值；根据前一棵 h.clust 树的高度设置纵轴范围。
     cex.lab=1.0,  # 轴标题按默认大小显示。
     cex.axis=1.0,  # 坐标轴刻度文字按默认大小显示。
     cex = 1,  # 样本标签按默认大小显示。
     cex.main=1  # 图标题按默认大小显示。
)  # 结束本次函数调用。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

# 去除可能的异常样本

# ---- 16. 切树筛选样本 ----
# 目的：按示例树高和最小簇大小生成样本簇，保留指定簇。
abline(h = 15, col = "red")  # abline() 在当前图的高度 15 处添加红色参考线。
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)  # cutreeStatic() 在高度 15 切树，最小簇为 10；当前 8 个样本可能全被排除。
table(clust)  # table() 统计每个簇中的样本数。
keepSamples = (clust==1)  # 生成逻辑索引，仅保留编号为 1 的样本簇。
protein2 = protein1[keepSamples, ]  # 按样本簇筛选表达表，保留全部蛋白列。
nprotein = ncol(protein2)  # ncol() 记录筛选后蛋白数，即表达表列数。
nSamples = nrow(protein2)  # nrow() 记录筛选后样本数，即表达表行数。

# 箱线图####

# ---- 17. 箱线图长表整理 ----
# 目的：将表达宽表展开成样本—表达值记录，并转换为 log2 尺度。
library(ggplot2)  # library() 加载 ggplot2，用于按图层绘图。
library(dplyr)  # library() 加载 dplyr，用于整理、筛选和分组数据。

# 转换长数据（接续之前的转换结果）
protein_long <- read.csv("4.箱线图/箱线图.csv", row.names = 1) %>%  # read.csv() 读取箱线图表达表，首列作行名，再传给下一步。
  tibble::rownames_to_column("protein") %>%  # tibble::rownames_to_column() 将蛋白行名转成名为 protein 的普通列。
  tidyr::pivot_longer(-protein, names_to = "sample", values_to = "expression")  # tidyr::pivot_longer() 宽表转长表，样本名存入 sample，表达值存入 expression。

protein_long_trans <- protein_long %>%  # 将长表传入后续转换，保存为新的数据对象。
  mutate(  # mutate() 开始新增对数表达列。
    expression_log2 = log2(expression + 1)  # log2() 计算 log2(expression+1)，加 1 以兼容零值。
  )  # 结束本次函数调用。
protein_long_trans = protein_long_trans[,-3]  # 删除原始表达列，保留对数表达列。

# 绘制箱线图

# ---- 18. 样本表达箱线图 ----
# 目的：对照各样本表达分布的位置、离散程度和离群值。
png("boxplot.png",width = 7,height = 6,units = "in",res = 300)  # png() 打开 7×6 英寸、300 dpi 的 PNG 图形设备。
pdf("boxplot.pdf",width = 7,height = 6)  # pdf() 打开 7×6 英寸的 PDF 图形设备。
ggplot(protein_long_trans, aes(x = sample, y = expression_log2)) +  # ggplot() 以样本为横轴、log2 表达为纵轴绘图；aes() 将数据变量映射到坐标或外观。
  geom_boxplot(  # geom_boxplot() 添加箱线图，后续参数控制颜色与样式。
    aes(  # aes() 开始设置箱线图的颜色映射。
      fill = substr(sample,1,1),  # substr() 用样本名首字母 H/L 映射填充颜色。
      color = substr(sample,1,1)  # substr() 用样本名首字母 H/L 映射边框颜色。
    ),  # 结束本次函数调用。
    # 关键参数调整
    coef = 1.5,  # 箱线须按 1.5 倍 IQR 判定，不是显示真实极值。
    fatten = 0.4,  # 将中位线相对粗细设为 0.4。
    size = 0.3,  # 将箱体边框线宽设为 0.3。
    outlier.shape = 19  # 用实心圆显示离群值，19 并非隐藏离群点。
  ) +  # 结束本次设置，继续添加图层。
# 补充箱线须端横线，统计规则仍沿用箱线图定义，不是改为真实极值。
  stat_boxplot(  # stat_boxplot() 添加箱线图须端的横线。
    aes(color = substr(sample,1,1)),  # aes() 须端横线按样本名首字母分组着色；substr() 按起止位置截取字符。
    geom = "errorbar",  # 将须端图层画成 errorbar 样式。
    width = 0.2,  # 设置须端横线宽度为 0.2。
    size = 0.5  # 设置须端横线线宽为 0.5。
  ) +  # 结束本次设置，继续添加图层。
  scale_fill_manual(  # scale_fill_manual() 手动设置箱体填充色及图例。
    name = "Group",  # 将图例标题设为 Group。
    values = c("H" = "#F6AB54", "L" = "#ADD8E6"),  # c() 将多个值合成向量；将 H 组填充为橙色，L 组填充为浅蓝色。
    labels = c("H" = "High", "L" = "Low")  # c() 将多个值合成向量；将图例中的 H/L 显示为 High/Low。
  ) +  # 结束本次设置，继续添加图层。
  scale_color_manual(  # scale_color_manual() 手动设置箱体边框颜色。
    name = "Group",  # 将图例标题设为 Group。
    values = c("H" = "#CE9759", "L" = "#549BB3"),  # c() 将多个值合成向量；将 H/L 边框分别设为深橙色和深蓝色。
    guide = "none"  # 隐藏边框颜色图例，避免重复。
  ) +  # 结束本次设置，继续添加图层。
  labs(title = "Protein Expression Distribution by Sample",  # labs() 设置样本表达分布图的标题。
       x = "Sample ID",  # 将横轴标题设为 Sample ID。
       y = "Expression Level (log2)") +  # 将纵轴标题设为 Expression Level (log2)。
  theme_minimal() +  # theme_minimal() 使用简洁主题。
  theme(  # theme() 继续设置图形主题。
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # element_text() 横轴标签旋转 45 度，右对齐并加粗。
    panel.grid.major.x = element_blank(),  # element_blank() 去掉竖向主要网格线。
    panel.border = element_rect(color = "#E0E0E0", fill = NA)  # element_rect() 添加浅灰色面板边框，不填充背景。
  )  # 结束本次函数调用。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

# 利用蛋白表达数据进行差异表达分析--------------------------------------------------

# ---- 19. 差异分析输入 ----
# 目的：读取表达表并执行当前脚本的缺失处理及 log2 转换。
exp = read.csv("5.差异分析/H_a-vs-L_a-exp.csv", row.names = 1)  # read.csv() 读取差异分析表达表，首列作为蛋白行名。

exp[is.na(exp)] <- 0  # is.na() 检查缺失值；将缺失值替换为 0；这不是按组均值填补。
exp = log2(exp+1)  # log2() 对表达矩阵计算 log2(x+1)，压缩数值范围。

# 构建分组矩阵

# ---- 20. 建立样本设计矩阵 ----
# 目的：用两组样本标签编码线性模型，确保分组与表达列顺序对应。
group = c(rep("H_a",4),rep("L_a",4))  # c() 将多个值合成向量；rep() 重复指定内容；构造前四个 H_a、后四个 L_a 标签，须与样本顺序一致。
dim(exp)  # dim() 查看数据行数和列数，核对输入规模。
design = model.matrix(~0 + factor(group))  # model.matrix() 建立无截距设计矩阵，每组各占一列；factor() 将类别编码为因子。
#model.matrix()根据分组信息生成一个模拟矩阵，线性拟合函数需要用模型矩阵
head(design)  # head() 查看设计矩阵列名，确认后续组名与列顺序对应。
colnames(design) = c("H_a", "L_a")  # colnames() 读取或设置列名；c() 将多个值合成向量；按实际组别顺序将两列命名为 H_a 和 L_a。
rownames(design) = colnames(exp)  # rownames() 读取或设置行名；colnames() 读取或设置列名；将设计矩阵行名设为表达表样本名，分组顺序仍须预先对齐。
head(design)  # head() 查看数据前六行，核对结构和内容。

# 进行差异表达分析

# ---- 21. 拟合 H_a−L_a 对比 ----
# 目的：对每个蛋白估计 H_a 相对 L_a 的表达变化。
fit = lmFit(exp, design)  # lmFit() 对每个蛋白拟合分组线性模型。
#design是一个设计矩阵，用于指定模型中每个样本的实验条件。它通常是一个因子向量或指示变量矩阵
cont.matrix = makeContrasts(contrasts="H_a-L_a", levels=design)  # makeContrasts() 设置 H_a−L_a 对比，正 logFC 表示 H_a 组更高。
fit2 = contrasts.fit(fit, cont.matrix)  # contrasts.fit() 按指定对比计算系数及标准误。

# 计算统计量和重要蛋白表

# ---- 22. 差异检验与结果提取 ----
# 目的：用经验贝叶斯稳定方差，并整理效应与 P 值。
fit2 = eBayes(fit2)  # eBayes() 用经验贝叶斯方法收缩方差估计，提高检验稳定性。
DEG_data1 = topTable(fit2, adjust="fdr", sort.by="B", number=nrow(exp))  # topTable() 提取全部蛋白结果，BH 校正 P 值，并按 B 统计量排序；nrow() 返回行数。
#提取差异分析结果：sort.by="B"按照B统计量（也称为lods或log-odds）对蛋白进行排序

# DEG 数据注释，# 筛选标准为p-value<0.05 & FC≥1.5或FC≤0.66

# ---- 23. 差异蛋白分类 ----
# 目的：按原始 P 值及 FC≥1.5 或 FC≤0.66 标记上下调。
DEG_data1$type = ifelse((DEG_data1$logFC >= log2(1.5) & DEG_data1$P.Value < 0.05),"up",  # ifelse() FC≥1.5 且原始 P<0.05 标为 up；这里未用校正 P 值筛选；log2() 取以 2 为底的对数。
                        ifelse((DEG_data1$logFC<= log2(0.66) & DEG_data1$P.Value < 0.05) ,"down","none"))  # ifelse() FC≤0.66 且原始 P<0.05 标为 down，其余为 none；log2() 取以 2 为底的对数。

DEG_data1$type1 = ifelse((DEG_data1$logFC >= log2(1.5) & DEG_data1$P.Value < 0.05),"up",  # ifelse() FC≥1.5 且原始 P<0.05 标为 up；这里未用校正 P 值筛选；log2() 取以 2 为底的对数。
                         ifelse((DEG_data1$logFC<= log2(0.66) & DEG_data1$P.Value < 0.05) ,"down","none"))  # ifelse() FC≤0.66 且原始 P<0.05 标为 down，其余为 none；log2() 取以 2 为底的对数。


# ---- 24. 设置分类标签并保存 ----
# 目的：统计类别数量，整理图例顺序，保存完整差异表。
lab_title = c(paste0("up:",length(which(DEG_data1$type=="up"))),  # c() 将多个值合成向量；paste0() 不加分隔符地拼接文字；length() 返回元素个数；which() 返回条件为 TRUE 的位置；统计上调蛋白数，并拼接 up 标签。
              paste0("no change:",length(which(DEG_data1$type=="none"))),  # paste0() 统计未达筛选阈值的蛋白数，并拼接 no change 标签；length() 返回元素个数；which() 返回条件为 TRUE 的位置。
              paste0("down:",length(which(DEG_data1$type=="down"))))  # paste0() 统计下调蛋白数，并拼接 down 标签；length() 返回元素个数；which() 返回条件为 TRUE 的位置。

DEG_data1$type = factor(DEG_data1$type,  # factor() 将 type 转为因子，用后续参数设置顺序和标签。
                        levels = c('up','none','down'),  # c() 将多个值合成向量；固定图例顺序为 up、none、down。
                        labels = lab_title)  # 用带蛋白数量的文字替换因子标签。

write.csv(DEG_data1,"protein.csv")  # write.csv() 保存完整蛋白差异结果及分类信息。

# 差异蛋白统计
# 加载必要包

# ---- 25. 汇总上下调数量 ----
# 目的：准备统计柱图；当前示例数量 57 和 15 为手工填写。
library(dplyr)  # library() 加载 dplyr，用于整理、筛选和分组数据。
library(ggplot2)  # library() 加载 ggplot2，用于按图层绘图。

diff_data=DEG_data1  # 复制差异分析结果，供后续整理使用。
#diff_data <- read.csv("protein1.csv",header = T)  # 备用代码，未执行；read.csv() 读取 CSV 表格。

# 统计上下调总数
summary_data <- data.frame(  # data.frame() 创建数据框；创建上下调蛋白数量汇总表。
  type1 = c("up","down"),  # c() 将多个值合成向量；定义汇总类别为 up 和 down。
  Total = c(57,15)  # c() 将多个值合成向量；手工填写上调 57、下调 15，未从当前结果自动统计。
)  # 结束本次函数调用。

# 修改因子顺序
summary_data <- summary_data %>%  # 将汇总表传给下一步整理。
  mutate(type1 = factor(type1, levels = c("up", "down")))  # mutate() 将类别设为因子，固定顺序为 up、down；factor() 将类别编码为因子；c() 将多个值合成向量。

# 基础绘图

# ---- 26. 单向数量柱形图 ----
# 目的：用柱高比较上下调蛋白数量，并在柱顶标注数值。
p <- ggplot(summary_data, aes(x = type1, y = Total, fill = type1)) +  # ggplot() 以调控方向为横轴、蛋白数量为纵轴绘制柱形图；aes() 将数据变量映射到坐标或外观。
  geom_col(width = 0.6, show.legend = FALSE) +  # geom_col() 柱宽为 0.6，隐藏填充图例。
  geom_text(aes(label = Total), vjust = -0.5, size = 5) +  # geom_text() 在柱顶添加数量标签，向上偏移 0.5，字号为 5；aes() 将数据变量映射到坐标或外观。
  scale_fill_manual(values = c("down" = "#2c7bb6", "up" = "#d7191c")) +  # scale_fill_manual() 将下调柱设为蓝色，上调柱设为红色；c() 将多个值合成向量。
  labs(title = "",  # labs() 图标题设为空。
       x = "Regulation Direction",  # 将横轴标题设为 Regulation Direction。
       y = "Number Of Proteins") +  # 将纵轴标题设为 Number Of Proteins。
  theme_minimal(base_size = 14) +  # theme_minimal() 使用简洁主题，基础字号为 14。
  theme(  # theme() 继续设置图形主题。
    panel.grid.major.x = element_blank(),  # element_blank() 去掉竖向主要网格线。
    axis.title = element_text(face = "bold"),  # element_text() 将坐标轴标题加粗。
    plot.title = element_text(hjust = 0.5, size = 16)  # element_text() 图标题居中，字号为 16。
  )  # 结束本次函数调用。

# 显示图形
print(p)  # print() 显示图形，并输出到当前图形设备。

# 数据预处理（添加坐标映射列）

# ---- 27. 双向数量柱形图 ----
# 目的：将上调画在零线上方、下调画在下方，突出变化方向。
summary_data <- summary_data %>%  # 将汇总表传入下一步，准备生成双向柱图坐标。
  mutate(Value = ifelse(type1 == "up", Total, -Total))  # mutate() 上调数量取正值、下调数量取负值，存入 Value；ifelse() 按条件逐元素选择返回值。

# 生成双向柱形图
png("hist.png",width = 7,height = 6,units = "in",res = 300)  # png() 打开 7×6 英寸、300 dpi 的 PNG 图形设备。
pdf("hist.pdf",width = 7,height = 6)  # pdf() 打开 7×6 英寸的 PDF 图形设备。
ggplot(summary_data, aes(x = type1, y = Value, fill = type1)) +  # ggplot() 以 Value 为柱高，按上下调方向着色；aes() 将数据变量映射到坐标或外观。
  geom_col(width = 0.4) +  # geom_col() 设置双向柱形图的柱宽为 0.4。
  scale_y_continuous(  # scale_y_continuous() 开始设置纵轴刻度、标签和留白。
    breaks = seq(-20, 60, by = 20),  # seq() 从 −20 到 60，每隔 20 设置一个纵轴刻度。
    labels = abs,  # abs() 取绝对值；刻度标签取绝对值，使下调柱也显示正数数量。
    name = "Number Of Proteins",  # 将纵轴标题设为 Number Of Proteins。
    expand = expansion(mult = c(0.1, 0.1))  # expansion() 纵轴两端各增加 10% 留白；c() 将多个值合成向量。
  ) +  # 结束本次设置，继续添加图层。
  geom_text(  # geom_text() 添加柱形图的数量标签。
    aes(label = Total,  # aes() 标签内容使用实际数量 Total。
        y = ifelse(type1 == "up", Value + 5, Value - 5)),  # ifelse() 上调标签上移 5，下调标签下移 5。
    color = "black",  # 将数量标签设为黑色。
    size = 4,  # 数量标签字号设为 4。
    fontface = "bold"  # 将数量标签加粗。
  ) +  # 结束本次设置，继续添加图层。
  scale_fill_manual(values = c("up" = "#e63946", "down" = "#1d3557")) +  # scale_fill_manual() 手动设置上调为红色，下调为蓝色；c() 将多个值合成向量。
  labs(x = "Regulation Direction",fill = "type") +  # labs() 设置横轴标题和填充图例标题。
  theme_minimal(base_size = 12) +  # theme_minimal() 使用简洁主题，基础字号为 12。
  theme(  # theme() 继续设置图形主题。
    panel.grid = element_blank(),  # element_blank() 去掉绘图区的全部网格线。
    legend.position = "right"  # 将图例放在右侧。
  ) +  # 结束本次设置，继续添加图层。
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +  # geom_hline() 在 y=0 处添加黑色实线作为基准线。
  geom_vline(xintercept = 0.4, color = "black", linetype = "solid", linewidth = 1)  # geom_vline() 在 x=0.4 处添加黑色竖线。
# solid 实线;dotted 点线;dotdash 点划线
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

# 差异蛋白火山图----------------------------------------------------------------

# ---- 28. 火山图结果输入 ----
# 目的：读取差异表，将效应大小与原始 P 值放在同一图中。
DEG_data1 = read.csv("5.差异分析/protein.csv", row.names = 1, header = T)  # read.csv() 读取已保存的蛋白差异结果，首列作为行名。

p_DEG = ggplot(DEG_data1, aes(logFC,-log10(P.Value),color = type))  # ggplot() 建立火山图对象，颜色按 type 分类映射；aes() 将数据变量映射到坐标或外观；log10() 取以 10 为底的对数。

lab_title = c(paste0("up:",length(which(DEG_data1$type1=="up"))),  # c() 将多个值合成向量；paste0() 不加分隔符地拼接文字；length() 返回元素个数；which() 返回条件为 TRUE 的位置；统计上调蛋白数，并拼接 up 标签。
              paste0("no change:",length(which(DEG_data1$type1=="none"))),  # paste0() 统计未达筛选阈值的蛋白数，并拼接 no change 标签；length() 返回元素个数；which() 返回条件为 TRUE 的位置。
              paste0("down:",length(which(DEG_data1$type1=="down"))))  # paste0() 统计下调蛋白数，并拼接 down 标签；length() 返回元素个数；which() 返回条件为 TRUE 的位置。

#火山图

# ---- 29. 基础火山图 ----
# 目的：用分类颜色和阈值线展示蛋白差异方向及检验证据。
p<- ggplot(DEG_data1,aes(x=logFC,y=-log10(P.Value)))+  # ggplot() 火山图横轴为 logFC，纵轴为原始 P 值的 −log10；aes() 将数据变量映射到坐标或外观；log10() 取以 10 为底的对数。
  scale_x_continuous(limits=c(-3,3), breaks=seq(-3,3,1))+  # scale_x_continuous() 横轴范围限定为 −3 至 3，每隔 1 设置刻度；c() 将多个值合成向量；seq() 生成等间距数列。
  geom_point(aes(color=type),size =1.5)+  # geom_point() 按 type 给蛋白散点着色，点大小为 1.5；aes() 将数据变量映射到坐标或外观。
  scale_color_manual(values=c("blue","gray","red"),labels = c("Down", "No Change", "Up"))+  # scale_color_manual() 按类别顺序分配蓝、灰、红色，并设置图例文字；c() 将多个值合成向量。
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8,col="#4C5B61")+  # geom_hline() 在原始 P=0.05 对应高度添加水平虚线；log10() 取以 10 为底的对数。
  geom_vline(xintercept = c(log2(1.5),log2(0.66)),lty=4,lwd=0.6,alpha=0.8,col="#4C5B61")+  # geom_vline() 添加 FC=1.5 和 0.66 对应的 log2FC 参考线；c() 将多个值合成向量；log2() 取以 2 为底的对数。
  annotate("text",x=2,y=5,label=lab_title[1],size=5,col="black") +  # annotate() 在 (2,5) 标注上调蛋白数量。
  annotate("text",x=-2,y=5,label=lab_title[3],size=5,col="black") +  # annotate() 在 (−2,5) 标注下调蛋白数量。
  theme_bw()+  # theme_bw() 使用白底带边框的主题。
  theme(panel.grid=element_blank(),axis.text = element_text(size = 12),axis.title = element_text(size = 12),legend.text = element_text(size = 12))  # theme() 去掉网格，轴刻度、轴标题及图例文字字号均设为 12；element_blank() 隐藏相应主题元素；element_text() 设置主题文字样式。
png("vloc_protein.png",width = 7,height = 6,units = "in",res = 300)  # png() 打开 7×6 英寸、300 dpi 的 PNG 图形设备。
pdf("vloc_protein.pdf",width = 7,height = 6)  # pdf() 打开 7×6 英寸的 PDF 图形设备。
p  # 显示图形对象 p；在脚本中运行时可显式 print()。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。


# ---- 30. 突出指定蛋白 ----
# 目的：在火山图上将五个指定蛋白用不同点样式突出显示。
proteinList0 <- c("RBM4B","NCAM1","PSAP","KLC1","CLIC1")  # c() 将多个值合成向量；指定待突出显示的五个蛋白名称。
proteinList <- DEG_data1[proteinList0,]  # 按蛋白行名提取指定的差异结果。

p<- ggplot(DEG_data1,aes(x=logFC,y=-log10(P.Value)))+  # ggplot() 火山图横轴为 logFC，纵轴为原始 P 值的 −log10；aes() 将数据变量映射到坐标或外观；log10() 取以 10 为底的对数。
  scale_x_continuous(limits=c(-5,5), breaks=seq(-5,5,1))+  # scale_x_continuous() 横轴范围限定为 −5 至 5，每隔 1 设置刻度；c() 将多个值合成向量；seq() 生成等间距数列。
  geom_point(aes(color=type),size=1.5)+  # geom_point() 按 type 给蛋白散点着色，点大小为 1.5；aes() 将数据变量映射到坐标或外观。
  scale_color_manual(values=c("blue","gray","red"),labels = c("Down", "No Change", "Up"))+  # scale_color_manual() 按类别顺序分配蓝、灰、红色，并设置图例文字；c() 将多个值合成向量。
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8,col="#4C5B61")+  # geom_hline() 在原始 P=0.05 对应高度添加水平虚线；log10() 取以 10 为底的对数。
  geom_vline(xintercept = c(log2(1.2),log2(0.83)),lty=4,lwd=0.6,alpha=0.8,col="#4C5B61")+  # geom_vline() 添加 FC=1.2、0.83 参考线，与前面的分类阈值不同；c() 将多个值合成向量；log2() 取以 2 为底的对数。
  geom_point(data=proteinList,aes(x = logFC, y = -log10(P.Value)),colour="yellow",size=2)+  # geom_point() 将选定的五个蛋白用黄色大点突出显示；aes() 将数据变量映射到坐标或外观；log10() 取以 10 为底的对数。
  annotate("text",x=3,y=11,label=lab_title[1],size=5,col="black") +  # annotate() 在 (3,11) 标注上调蛋白数量。
  annotate("text",x=-3,y=11,label=lab_title[3],size=5,col="black") +  # annotate() 在 (−3,11) 标注下调蛋白数量。
  theme_bw()+  # theme_bw() 使用白底带边框的主题。
  theme(panel.grid=element_blank(),axis.text = element_text(size = 12),axis.title = element_text(size = 12),legend.text = element_text(size = 12))  # theme() 去掉网格，轴刻度、轴标题及图例文字字号均设为 12；element_blank() 隐藏相应主题元素；element_text() 设置主题文字样式。

pdf("vloc.pdf",width = 7,height = 6)  # pdf() 打开 7×6 英寸的 PDF 图形设备。
png("vloc.png",width = 6,height = 5,units = "in",res = 300)  # png() 打开 6×5 英寸、300 dpi 的 PNG 图形设备。
p  # 显示图形对象 p；在脚本中运行时可显式 print()。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。
write.csv(df,"df.csv")  # write.csv() 导出 df；本脚本此前未定义该对象，运行前需核对。
#标记出5个蛋白的label

# ---- 31. 火山图名称标注 ----
# 目的：匹配五个蛋白名称并添加自动避让的文字标签。
proteinList1 <- DEG_data1[rownames(DEG_data1) %in% proteinList0,]  # rownames() 读取或设置行名；按名称匹配五个蛋白，仅保留实际存在的记录。
proteinList1 <- subset(proteinList1, select = -type)  # subset() 删除 type 列，避免继承原图的颜色分类映射。
proteinList1$label <- rownames(proteinList1)  # rownames() 读取或设置行名；把蛋白行名存入 label 列，供添加文字标签。

library(ggrepel)  # library() 加载 ggrepel，用于添加自动避让的文字标签。
pdf("vloc_protein_label.pdf",width = 6,height = 5)  # pdf() 打开 6×5 英寸的 PDF 图形设备。
png("vloc_protein_label.png",width = 6,height = 5,units = "in",res = 300)  # png() 打开 6×5 英寸、300 dpi 的 PNG 图形设备。
p + geom_label_repel(data = proteinList1,  # geom_label_repel() 添加自动避让的标签框；在原火山图上叠加自动避让的蛋白标签。
                     aes(x = logFC, y = -log10(P.Value), label = label),  # aes() 使用 logFC、−log10(P) 定位，标签内容取 label 列；log10() 取以 10 为底的对数。
                     size = 3,color="black",  # 将标签字号设为 3，文字设为黑色。
                     box.padding = unit(0.4, "lines"),  # unit() 给长度附加单位；设置标签框周围留白为 0.4 行高。
                     segment.color = "black",  # 将标签与散点之间的连线设为黑色。
                     segment.size = 0.4,  # 将标签连线线宽设为 0.4。
                     max.overlaps = 20  # 标签允许的重叠数量上限设为 20。
)  # 结束本次函数调用。

dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

# 聚类分析热图------------------------------------------------------------------
#BiocManager::install("ggvenn")  # 备用代码，未执行；BiocManager::install() 安装或更新软件包。

# ---- 32. 全部差异蛋白表达 ----
# 目的：提取上下调差异蛋白及其样本表达，为候选热图准备数据。
library(ComplexHeatmap)  # library() 加载 ComplexHeatmap，用于绘制带注释的热图。
library(pheatmap)  # library() 加载 pheatmap，用于绘制聚类热图。

# 提取全部差异蛋白作为候选热图输入；不是展示所有被检测蛋白。
down_rows <- DEG_data1[DEG_data1$type1 == "down", ]  # 提取 type1 为 down 的差异蛋白。
up_rows <- DEG_data1[DEG_data1$type1 == "up",]  # 提取 type1 为 up 的差异蛋白。
up_down <- rbind(down_rows,up_rows)  # rbind() 按行合并下调与上调蛋白结果。
write.csv(up_down,"up_down_FC.csv")  # write.csv() 保存全部上下调差异蛋白结果。
exp_up_down_protein <- exp[rownames(up_down),]  # rownames() 读取或设置行名；按差异蛋白行名提取其样本表达值。

#按照蛋白差异倍数进行排序

# ---- 33. 按效应挑选候选 ----
# 目的：选择上调 20 个和下调 15 个蛋白，形成变化幅度较大的子集。
protein_FC <- DEG_data1[order(DEG_data1$logFC),]  # order() 返回排序后的元素位置；将蛋白按 logFC 从小到大排序。
down_20protein <- protein_FC[grep("down",protein_FC$type1)[1:15],]  # grep() 返回匹配字符串的位置；取 logFC 最小的 15 个下调蛋白；对象名虽写 20，实际取 15。
up_20protein <- tail(protein_FC[grep("up",protein_FC$type1),], n =20)  # tail() 取 logFC 最大的 20 个上调蛋白；grep() 返回匹配字符串的位置。
updown_20protein <- rbind(up_20protein,down_20protein)  # rbind() 合并所选上调 20 个、下调 15 个蛋白。
exp_updown_20protein <- exp[rownames(updown_20protein),]  # rownames() 读取或设置行名；按所选蛋白名称提取表达矩阵。

# 加入某一蛋白

# ---- 34. 加入指定候选 ----
# 目的：将 SCN4A 和 NCAM1 加入另存表达表；后续热图仍使用原候选表。
selected_protein <- exp[c("SCN4A","NCAM1"), ]  # c() 将多个值合成向量；额外提取 SCN4A 和 NCAM1 的表达值。
selected_all_protein <- rbind(selected_protein,exp_updown_20protein)  # rbind() 将额外选定蛋白与上下调候选表达表按行拼接。
write.csv(selected_all_protein,"protein_total_FC.csv")  # write.csv() 保存加入指定蛋白后的表达表。

#BiocManager::install("ComplexHeatmap")  # 备用代码，未执行；BiocManager::install() 安装或更新软件包。

# ---- 35. ComplexHeatmap 热图 ----
# 目的：对当前候选矩阵做行内标准化，展示组内及组间相对表达模式。
library(ComplexHeatmap)  # library() 加载 ComplexHeatmap，用于绘制带注释的热图。
C <- exp_updown_20protein  # 使用原上下调候选表作热图，未包含刚额外加入的蛋白。
C <- as.matrix(C)  # as.matrix() 转为矩阵；将热图输入转为数值矩阵。
samples <- c(rep("H_a",4),rep("L_a",4))  # c() 将多个值合成向量；rep() 重复指定内容；构造前四个 H_a、后四个 L_a 标签，须与样本顺序一致。
for (i in 1:nrow(C)) C[i, ] <- scale(log(unlist(C[i, ] + 1), 2))  # nrow() 返回行数；scale() 默认按列减均值并除以标准差；log() 取指定底数的对数；unlist() 将列表展开为向量；按蛋白计算 log2(x+1) 后做 z-score；exp 已取过对数，会重复变换。
C[is.na(C)]=0  # is.na() 检查缺失值；将标准化后的 NA 替换为 0。
#C <- apply(C, 1, scale) #标准化方式二  # 备用代码，未执行；apply() 沿指定行或列执行函数。
D <- Heatmap(C,  # Heatmap() 使用矩阵 C 创建 ComplexHeatmap 热图对象。
             name = "expression",  # 将热图图例标题设为 expression。
             col = colorRampPalette(c("navy","white","firebrick3"))(100),  # colorRampPalette() 创建渐变配色函数；c() 将多个值合成向量；用深蓝—白—深红的 100 色渐变表示数值高低。
             show_row_names = T,  # 显示蛋白行名，TRUE 并非隐藏行名。
             show_column_names = F,  # 隐藏样本列名。
             column_names_gp = gpar(fontsize = 8),  # gpar() 设置 grid 图形的字体等属性；样本列名字号设为 8；隐藏时不显示。
             column_names_rot = 90,  # 将样本列名旋转 90 度。
             #cluster_columns = F, 
             row_names_gp = gpar(fontsize = 8),  # gpar() 设置 grid 图形的字体等属性；将蛋白行名字号设为 8。
             top_annotation = HeatmapAnnotation(Group = samples,  # HeatmapAnnotation() 在热图顶部添加样本分组注释。
                                                simple_anno_size = unit(2, 'mm'),  # unit() 给长度附加单位；将分组注释条厚度设为 2 mm。
                                                col = list(Group = c('L_a' = '#00DAE0',  # list() 将对象组成列表；c() 将多个值合成向量；在分组色表中指定 L_a 为青色。
                                                                     'H_a' = '#FF9289')),  # 指定 H_a 为粉色，完成分组色表。
                                                show_annotation_name = FALSE))  # 隐藏顶部注释名称，仅保留色条及图例。


# ---- 36. 保存第一版热图 ----
# 目的：将热图对象输出到当前图形设备。
pdf("heatmap_FC.pdf",width = 7,height = 6)  # pdf() 打开 7×6 英寸的 PDF 图形设备。
png("heatmap_FC.png",width = 7,height = 6,units = "in",res = 300)  # png() 打开 7×6 英寸、300 dpi 的 PNG 图形设备。
D  # 显示图形对象 D；在脚本中运行时可显式 print()。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

# 聚类热图可视化----------

# ---- 37. 重建聚类热图输入 ----
# 目的：重新提取按效应挑选的上下调候选，准备第二种热图展示。
library(ComplexHeatmap)  # library() 加载 ComplexHeatmap，用于绘制带注释的热图。
library(pheatmap)  # library() 加载 pheatmap，用于绘制聚类热图。

#按照蛋白差异倍数进行排序
protein_FC <- DEG_data1[order(DEG_data1$logFC),]  # order() 返回排序后的元素位置；将蛋白按 logFC 从小到大排序。
down_20protein <- protein_FC[grep("down",protein_FC$type1)[1:15],]  # grep() 返回匹配字符串的位置；取 logFC 最小的 15 个下调蛋白；对象名虽写 20，实际取 15。
up_20protein <- tail(protein_FC[grep("up",protein_FC$type1),], n =20)  # tail() 取 logFC 最大的 20 个上调蛋白；grep() 返回匹配字符串的位置。
updown_20protein <- rbind(up_20protein,down_20protein)  # rbind() 合并所选上调 20 个、下调 15 个蛋白。

#提取表达量矩阵
exp_updown_20protein <- exp[rownames(updown_20protein),]  # rownames() 读取或设置行名；按所选蛋白名称提取表达矩阵。

library(pheatmap)  # library() 加载 pheatmap，用于绘制聚类热图。
library(RColorBrewer)  # library() 加载 RColorBrewer，用于调用预设调色板。
library(ComplexHeatmap)  # library() 加载 ComplexHeatmap，用于绘制带注释的热图。

# 数据标准化（按行Z-score）

# ---- 38. 按蛋白行标准化 ----
# 目的：将每个蛋白跨样本变换为 z-score，使不同蛋白的相对模式可比较。
scaled_matrix <- t(scale(t(exp_updown_20protein)))  # t() 转置行列；scale() 默认按列减均值并除以标准差；按蛋白行做 z-score：减去行均值，再除以行标准差。

# 创建注释数据框

# ---- 39. 热图分组注释 ----
# 目的：建立与样本列对应的组别及颜色，辅助读取样本聚类结果。
annotation_df <- data.frame(  # data.frame() 创建数据框；创建热图样本注释数据框。
  Group = factor(c(rep("H_a",4),rep("L_a",4)))  # factor() 按四个 H_a、四个 L_a 建立分组因子；c() 将多个值合成向量；rep() 重复指定内容。
)  # 结束本次函数调用。

rownames(annotation_df) <- colnames(scaled_matrix)  # rownames() 读取或设置行名；colnames() 读取或设置列名；用表达矩阵样本名设置注释行名，保证身份对应。

# 定义注释颜色
annotation_colors <- list(  # list() 将对象组成列表；创建分组注释的颜色列表。
  Group = c(H_a="#66C2A5", L_a="#E1917A")  # c() 将多个值合成向量；将 H_a 注释设为青绿色，L_a 设为粉橙色。
  #Batch = c("1"="gray90", "2"="gray40")  # 备用代码，未执行；c() 将多个值合成向量。
)  # 结束本次函数调用。

# 绘制聚类热图

# ---- 40. 聚类热图样式 ----
# 目的：固定蛋白行顺序并聚类样本列，用注释色条展示组别。
library(conflicted)  # library() 加载 conflicted，用于指定同名函数优先使用的来源。
conflicts_prefer(ComplexHeatmap::pheatmap)  # conflicts_prefer() 指定优先使用 ComplexHeatmap 的 pheatmap 实现；ComplexHeatmap::pheatmap() 创建聚类热图。

D <- pheatmap(  # pheatmap() 创建聚类热图，后续逐项设置参数。
  mat = scaled_matrix,  # 使用按蛋白行标准化的表达矩阵。
  #color = colorRampPalette(rev(brewer.pal(n=11, name="RdYlBu")))(100),  # 备用代码，未执行；colorRampPalette() 创建渐变配色函数；rev() 反转向量顺序；brewer.pal() 从命名调色板提取颜色。
  col = colorRampPalette(c("navy","white","firebrick3"))(100),  # colorRampPalette() 创建渐变配色函数；c() 将多个值合成向量；用深蓝—白—深红的 100 色渐变表示数值高低。
  cluster_rows = F,  # 关闭蛋白行聚类，保留当前行顺序。
  cluster_cols = T,  # 对样本列进行聚类。
  clustering_method = "ward.D2",  # 聚类连接方法设为 ward.D2。
  show_rownames = T,  # 显示蛋白行名。
  show_colnames = F,  # 隐藏样本列名。
  annotation_col = annotation_df,  # 传入与样本列对应的分组注释表。
  annotation_colors = annotation_colors,  # 使用预先定义的分组注释颜色。
  border_color = "white",  # 将热图单元格边框设为白色。
  fontsize_col = 8,  # 样本列名字号设为 8。
  fontsize_row = 8,  # 蛋白行名字号设为 8。
  treeheight_row = 13,  # 设置行聚类树高度为 13；当前行聚类已关闭。
  cellwidth = 40,  # 将单元格宽度设为 40。
  cellheight = 10,  # 将单元格高度设为 10。
  main = "Clustered Heatmap of Gene Expression"  # 设置热图主标题。
)  # 结束本次函数调用。

# ---- 41. 保存第二版热图 ----
# 目的：输出带样本分组和列聚类的热图对象。
pdf("ComplexHeatmap_FC.pdf",width = 8,height = 8.5)  # pdf() 打开 8×8.5 英寸的 PDF 图形设备。
png("ComplexHeatmap_FC.png",width = 8,height = 8.5,units = "in",res = 300)  # png() 打开 8×8.5 英寸、300 dpi 的 PNG 图形设备。
D  # 显示图形对象 D；在脚本中运行时可显式 print()。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

#GO富集分析---------------------------------------------------------------------
# 加载富集分析所需要的R包 
#BiocManager::install("clusterProfiler")  # 备用代码，未执行；BiocManager::install() 安装或更新软件包。


# ---- 42. 加载富集工具 ----
# 目的：准备标识转换和功能富集；本例实际使用猪注释库。
library(clusterProfiler)  # library() 加载 clusterProfiler，用于进行标识转换和功能富集。
library(org.Hs.eg.db)  # library() 加载 org.Hs.eg.db，用于查询人类基因注释。
#BiocManager::install("org.Ss.eg.db")  # 备用代码，未执行；BiocManager::install() 安装或更新软件包。
library(org.Ss.eg.db)  # library() 加载 org.Ss.eg.db，用于查询猪基因注释。
library(ggplot2)  # library() 加载 ggplot2，用于按图层绘图。

# 读取数据：蛋白差异分析的结果 

# ---- 43. 富集成员标识转换 ----
# 目的：将候选蛋白的基因符号转换为猪 ENTREZID。
diffpro <- read.csv("6.富集分析/富集分析.csv",header = T)  # read.csv() 读取富集分析用的差异蛋白表。
rownames(diffpro) <- diffpro$X  # rownames() 读取或设置行名；将 X 列中的蛋白标识设为行名。
Protein <- bitr(rownames(diffpro),  # bitr() 开始转换蛋白对应基因的标识类型；rownames() 读取或设置行名。
                fromType = "SYMBOL",  # 输入标识为基因符号 SYMBOL。
                toType = c("ENTREZID"),  # c() 将多个值合成向量；将基因符号转换为 ENTREZID。
                OrgDb = org.Ss.eg.db)  # 使用猪注释数据库 org.Ss.eg.db。

# GO富集分析与结果存储

# ---- 44. GO 过度富集 ----
# 目的：比较候选基因与背景的 GO 成员比例，分析 BP、CC、MF 三类条目。
GO <- enrichGO(gene = Protein$ENTREZID,  # enrichGO() 用候选蛋白对应的 ENTREZID 进行 GO 过度富集分析。
               OrgDb = org.Ss.eg.db,  # 使用猪的 GO 注释数据库。
               keyType = "ENTREZID",  # 指定输入标识类型为 ENTREZID。
               ont = "all",  # 同时分析 BP、CC、MF 三类 GO 条目。
               pAdjustMethod = "BH",  # 使用 BH 方法校正多重检验 P 值。
               pvalueCutoff = 1,  # P 值阈值设为 1，保留结果供后续筛选。
               qvalueCutoff = 1,  # q 值阈值设为 1，不在此处收紧筛选。
               minGSSize = 5,  # 仅检验背景中至少含 5 个基因的条目，并非命中数至少为 5。
               maxGSSize = 2000,  # 仅检验背景中最多含 2000 个基因的条目。
               readable = TRUE)  # 将结果成员的 ENTREZID 转为可读基因符号。

# ---- 45. 保存 GO 结果 ----
# 目的：将条目 ID、描述、成员及统计量写入 CSV 和文本表。
write.csv(data.frame(ID=row.names(GO@result),GO@result),  # write.csv() 将 GO 结果转为数据框，并增加条目 ID 列；data.frame() 创建数据框；row.names() 读取或设置行名。
          file = "GO_result_Protein.csv")  # 指定 GO 结果 CSV 的保存文件名。

write.table(data.frame(ID=row.names(GO@result), GO@result),  # write.table() 将带条目 ID 的 GO 结果写为文本表格；data.frame() 创建数据框；row.names() 读取或设置行名。
            file = "GO_result_Protein.txt",  # 指定 GO 结果文本文件名。
            sep = "\t",  # 用制表符分隔各列。
            quote = FALSE,  # 输出字符内容时不额外加引号。
            row.names = FALSE)  # 不另行输出 R 行名。

# 结果可视化--柱状图

# ---- 46. GO 富集柱形图 ----
# 目的：按 GO 类别分面，对照各条目的命中数量。
pdf(file="go_barplot_pro.pdf",width = 6,height = 5)  # pdf() 打开 6×5 英寸的 PDF 图形设备。
png("go_barplot_pro.png",width = 6,height = 5,units = "in",res = 300)  # png() 打开 6×5 英寸、300 dpi 的 PNG 图形设备。
barplot(GO, drop = TRUE,  # barplot() 绘制 GO 富集柱形图，默认柱长表示命中数 Count。
        showCategory =10,  # 每个 GO 类别最多展示 10 个条目。
        label_format=100,  # 条目名称超过 100 字符时换行。
        split="ONTOLOGY") +  # 按 ONTOLOGY 区分 BP、CC、MF。
  facet_grid(ONTOLOGY~., scale='free') +  # facet_grid() 按 GO 类别纵向分面，各面板使用独立坐标尺度。
  theme(axis.text.x = element_text(size = 10),  # theme() 将横轴刻度字号设为 10；element_text() 设置主题文字样式。
        axis.text.y = element_text(size = 6),  # element_text() 将纵轴条目名字号设为 6。
        strip.text = element_text(size = 8))  # element_text() 将分面标题字号设为 8。
#panel.spacing = unit(1, "lines")) #调整面板之间的间距  # 备用代码，未执行；unit() 给长度附加单位。
#axis.title.x = element_text(size = 14),  # 调整x轴标题字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
#axis.title.y = element_text(size = 14),  # 调整y轴标题字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
#legend.title = element_text(size = 12),  # 调整图例标题字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
#legend.text = element_text(size = 10))   # 调整图例文本字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。
head(GO@result)  # head() 查看 GO 结果表的前六行。

# 结果可视化--点图

# ---- 47. GO 富集分面点图 ----
# 目的：用位置、点大小和颜色联合展示 GO 富集结果。
pdf(file="go_dotplot_pro.pdf",width = 6,height = 5)  # pdf() 打开 6×5 英寸的 PDF 图形设备。
png("go_dotplot_pro.png",width = 6,height = 5,units = "in",res = 300)  # png() 打开 6×5 英寸、300 dpi 的 PNG 图形设备。
# 方式一：分面
dotplot(GO,showCategory = 10,label_format=100,split="ONTOLOGY") +  # dotplot() 绘制 GO 点图，按类别拆分，每类最多展示 10 项。
  facet_grid(ONTOLOGY~., scale='free') +  # facet_grid() 按 GO 类别纵向分面，各面板使用独立坐标尺度。
  scale_size(range = c(2, 4)) +  # scale_size() 将点大小的显示范围设为 2 至 4；c() 将多个值合成向量。
  theme(  # theme() 继续设置图形主题。
    #text = element_text(size = 14),  # 修改整体字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
    #axis.title = element_text(size = 16),  # 修改轴标题字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
    axis.text.y = element_text(size = 6),  # element_text() 将纵轴条目名字号设为 6。
    axis.text.x = element_text(size = 8)  # element_text() 将横轴刻度字号设为 8。
    #legend.title = element_text(size = 14),  # 修改图例标题字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
    #legend.text = element_text(size = 12)   # 修改图例文本字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
  )  # 结束本次函数调用。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

# 方式二：不分面

# ---- 48. GO 不分面点图 ----
# 目的：在一个面板中按 GeneRatio 对照最多 15 个条目。
dotplot(  # dotplot() 开始绘制不分面的 GO 富集点图。
  GO,  # 传入 GO 富集结果对象。
  x = "GeneRatio",  # 横轴用 GeneRatio：命中数除以有效可注释输入成员数。
  color = "p.adjust",  # 用 BH 校正后 P 值映射颜色。
  showCategory = 15,  # 最多展示 15 个富集条目。
  size = NULL,  # 不手工指定点大小字段，采用默认映射。
  split = NULL,  # 不按其他字段拆分条目。
  font.size = 12,  # 图中文字基础字号设为 12。
  title = "",  # 图标题设为空。
  orderBy = "x",  # 按横轴数值对展示条目排序。
  label_format = 30)  # 条目名称超过 30 字符时换行。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

# KEGG富集分析与结果存储 -------------------------------------------------------

# ---- 49. KEGG 过度富集 ----
# 目的：用猪通路注释检验候选成员的 KEGG 通路富集。
KEGG <- enrichKEGG(gene = Protein$ENTREZID,  # enrichKEGG() 用猪基因 ENTREZID 进行 KEGG 过度富集分析。
                   organism = "ssc",  # 指定 KEGG 物种为猪 ssc。
                   keyType = "kegg",  # 使用 KEGG 的基因编号类型。
                   pAdjustMethod = "BH",  # 使用 BH 方法校正多重检验 P 值。
                   pvalueCutoff = 1,  # P 值阈值设为 1，保留结果供后续筛选。
                   qvalueCutoff = 1)  # q 值阈值设为 1，保留结果供后续筛选。


# ---- 50. KEGG 富集柱形图 ----
# 目的：以 GeneRatio 为柱长，对照条目命中比例和校正 P 值。
pdf(file="KEGG_barplot_pro.pdf",width = 6,height = 5)  # pdf() 打开 6×5 英寸的 PDF 图形设备。
png(file="KEGG_barplot_pro.png",width = 6,height = 5,res = 300,units = "in")  # png() 打开 6×5 英寸、300 dpi 的 PNG 图形设备。
barplot(KEGG,  # barplot() 绘制 KEGG 富集柱形图。
        x = "GeneRatio",  # 横轴用 GeneRatio：命中数除以有效可注释输入成员数。
        color = "p.adjust",  # 用 BH 校正后 P 值映射颜色。
        label_format=100,  # 条目名称超过 100 字符时换行。
        showCategory = 15,  # 最多展示 15 个富集条目。
        title = "KEGG_enrichment")+  # 将图标题设为 KEGG_enrichment。
  theme(  # theme() 继续设置图形主题。
    #text = element_text(size = 14),  # 修改整体字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
    #axis.title = element_text(size = 16),  # 修改轴标题字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
    axis.text.y = element_text(size = 8),  # element_text() 将纵轴条目名字号设为 8。
    axis.text.x = element_text(size = 8)  # element_text() 将横轴刻度字号设为 8。
    #legend.title = element_text(size = 14),  # 修改图例标题字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
    #legend.text = element_text(size = 12)   # 修改图例文本字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
  )  # 结束本次函数调用。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。


# ---- 51. KEGG 富集点图 ----
# 目的：用点图展示 KEGG 条目的命中数量与校正 P 值。
pdf(file="KEGG_dotplot_pro.pdf",width = 6,height = 5)  # pdf() 打开 6×5 英寸的 PDF 图形设备。
png("KEGG_dotplot_pro.png",width = 6,height = 5,units = "in",res = 300)  # png() 打开 6×5 英寸、300 dpi 的 PNG 图形设备。
dotplot(KEGG,color = "p.adjust",label_format=100,)+  # dotplot() 绘制 KEGG 点图，颜色取 p.adjust；末尾多余逗号可能影响运行。
  theme(  # theme() 继续设置图形主题。
    #text = element_text(size = 14),  # 修改整体字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
    #axis.title = element_text(size = 16),  # 修改轴标题字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
    axis.text.y = element_text(size = 8),  # element_text() 将纵轴条目名字号设为 8。
    axis.text.x = element_text(size = 8))+  # element_text() 将横轴刻度字号设为 8。
  #legend.title = element_text(size = 14),  # 修改图例标题字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
  #legend.text = element_text(size = 12)   # 修改图例文本字体大小  # 备用代码，未执行；element_text() 设置主题文字样式。
  scale_size(range = c(3,6))  # scale_size() 将点大小的显示范围设为 3 至 6；c() 将多个值合成向量。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

# 结果存储

# ---- 52. 可读 KEGG 结果与工作空间 ----
# 目的：转换结果成员名称、导出表格并保存当前 R 对象。
KEGG_results <- clusterProfiler::setReadable(KEGG,  # clusterProfiler::setReadable() 将 KEGG 结果中的成员标识转换为可读名称。
                                  OrgDb = org.Ss.eg.db,  # 用猪注释库匹配成员名称，不改变富集检验结果。
                                  keyType = "ENTREZID")  # 指定待转换的成员标识为 ENTREZID。
write.csv(data.frame(ID=rownames(KEGG_results@result),KEGG_results@result),  # write.csv() 将可读 KEGG 结果转为数据框，并增加条目 ID；data.frame() 创建数据框；rownames() 读取或设置行名。
          file = "KEGG_enrich_result_pro.csv")  # 指定 KEGG 结果 CSV 的保存文件名。

write.table(data.frame(ID=rownames(KEGG_results@result),KEGG_results@result),  # write.table() 将可读 KEGG 结果写为文本表格；data.frame() 创建数据框；rownames() 读取或设置行名。
            file = "KEGG_enrich_result.txt",  # 指定 KEGG 结果文本文件名。
            sep = "\t",  # 用制表符分隔各列。
            quote = FALSE,  # 输出字符内容时不额外加引号。
            row.names = FALSE)  # 不另行输出 R 行名。

save.image(file = "workspace.RData")  # save.image() 将当前工作空间中的全部对象保存为 RData。

#弦图---------------------------------------------------------------------------
# 加载所需要的包
#install_github("wencke/wencke.github.io")  # 备用代码，未执行；install_github() 从 GitHub 安装软件包。

# ---- 53. 弦图输入读取 ----
# 目的：读取通路富集表及蛋白差异表，准备连接成员和 logFC。
library(tidyverse)  # library() 加载 tidyverse，用于整理数据和绘图。
library(GOplot)  # library() 加载 GOplot，用于整理富集结果并绘制弦图。
library(RColorBrewer)  # library() 加载 RColorBrewer，用于调用预设调色板。
library(clusterProfiler)  # library() 加载 clusterProfiler，用于进行标识转换和功能富集。
KEGG <- read.csv("7.弦图/KEGG_enrich_result_pro.csv")  # read.csv() 读取弦图所用的 KEGG 富集结果表。
diffprotein <- read.csv("7.弦图/protein1.csv",header = T)  # read.csv() 读取蛋白差异结果，用于匹配 logFC。
#DEGlgFC <- read.csv("P_0.05_DEGs_genes.csv")  # 备用代码，未执行；read.csv() 读取 CSV 表格。
#DEGlgFC <- DEGlgFC[,c('gene','log2FoldChange')]  # 备用代码，未执行；c() 将多个值合成向量。


# ---- 54. 选择展示通路 ----
# 目的：查找感兴趣名称并按当前固定行号提取四个条目。
grep('Parkinson disease',KEGG$Description)  # grep() 查找 Parkinson disease 所在行号。
grep('Spinocerebellar ataxia',KEGG$Description)  # grep() 查找 Spinocerebellar ataxia 所在行号。
grep('Alzheimer disease',KEGG$Description)  # grep() 查找 Alzheimer disease 所在行号。
grep('Human T-cell leukemia virus 1 infection',KEGG$Description)  # grep() 查找指定病毒感染通路所在行号。
KEGG <- KEGG[c(3,5,8,12),c('Description','geneID')]  # c() 将多个值合成向量；固定取第 3、5、8、12 行的通路名称和成员列，须核对实际条目。
#KEGG_selected$Description <- sub(" - .*", "", KEGG_selected$Description)

head(KEGG$Description)  # head() 查看所选通路名称。


# ---- 55. 整理通路成员 ----
# 目的：拆分成员 ID，并按标识补充每个成员的 logFC。
KEGG_Clean = KEGG %>%  # 将通路表传入下一步，准备拆分成员列。
  separate_rows(geneID,sep = "/")  # separate_rows() 按 / 拆分 geneID，一行保留一个通路—成员对应关系。
LF <- data.frame(geneID = unique(KEGG_Clean$geneID))  # data.frame() 创建数据框；unique() 提取不重复值；提取不重复成员，创建用于匹配 logFC 的数据框。
LF$logFC <- diffprotein$logFC[match(LF$geneID,diffprotein$X)]  # match() 查找标识在参照表中的位置；按成员标识匹配差异表 X 列，加入 logFC；未匹配时为 NA。

KEGG_Clean <- KEGG_Clean %>%  # 将通路成员表传入下一步，准备连接差异倍数。
  left_join(LF,by="geneID")  # left_join() 按 geneID 左连接 logFC，保留全部通路成员。
colnames(KEGG_Clean) <- c('Term','Genes','logFC')  # colnames() 读取或设置列名；c() 将多个值合成向量；将三列统一命名为 Term、Genes、logFC。


# ---- 56. 构建弦图关联矩阵 ----
# 目的：每条通路保留当前前 10 个成员，转换为成员—通路关系矩阵。
library(conflicted)  # library() 加载 conflicted，用于指定同名函数优先使用的来源。
conflicts_prefer(clusterProfiler::slice)  # conflicts_prefer() 优先选 clusterProfiler::slice，需确认其适用于当前数据框；clusterProfiler::slice() 按位置选取数据行。
KEGG_Clean <- KEGG_Clean %>%  # 将含 logFC 的通路成员表传入分组步骤。
  group_by(Term) %>%  # group_by() 按通路 Term 分组。
  slice(1:10) %>%  # slice() 每个通路取当前顺序的前 10 行，未按 logFC 再排序。
  ungroup()  # ungroup() 取消分组，恢复普通表格。
dfPlot = chord_dat(data.frame(KEGG_Clean), process=unique(KEGG_Clean$Term))  # chord_dat() 构建成员—通路关联矩阵，供弦图使用；data.frame() 创建数据框；unique() 提取不重复值。

# 把log2FC添加进去

# ---- 57. 补充弦图数值 ----
# 目的：按矩阵行名匹配 logFC，保证成员颜色与差异方向对应。
new_col <- diffprotein$logFC[match(rownames(dfPlot),diffprotein$X)]  # match() 查找标识在参照表中的位置；rownames() 读取或设置行名；按弦图矩阵行名匹配蛋白 logFC。
dfPlot <- cbind(dfPlot, new_col)  # cbind() 将匹配的 logFC 按列拼接到弦图矩阵。
colnames(dfPlot)[ncol(dfPlot)] <- "logFC"  # colnames() 读取或设置列名；ncol() 返回列数；将新增的最后一列命名为 logFC。

#custom_colors <- c("#ce9fe4","#00a9e9","#d2e7d1","#d6d7f1","#ff7976") #可以设置自己喜欢的颜色  # 备用代码，未执行；c() 将多个值合成向量。

# ---- 58. 绘制功能弦图 ----
# 目的：展示所选通路与蛋白成员的归属关系，并按 logFC 着色。
library(ggplot2)  # library() 加载 ggplot2，用于按图层绘图。
library(ggrepel)  # library() 加载 ggrepel，用于添加自动避让的文字标签。
dfPlot1 = as.data.frame(dfPlot)  # as.data.frame() 转为数据框；将弦图输入转换为数据框。
GOChord(dfPlot1,  # GOChord() 绘制成员—通路弦图，连带表示成员属于相应通路。
        title="KEGG_enriched",  # 将弦图标题设为 KEGG_enriched。
        space = 0.05,  # 设置蛋白方格之间的间隙为 0.05。
        gene.order = 'logFC',  # 按 logFC 排列蛋白。
        gene.space = 0.25,  # 将蛋白标签与图形的间距设为 0.25。
        gene.size = 4,  # 蛋白标签字号设为 4。
        lfc.col=c('red','lightblue','darkblue'),  # c() 将多个值合成向量；设置 logFC 渐变配色，颜色高低方向以实际图例为准。
        ribbon.col=brewer.pal(ncol(dfPlot)-1, "Set3"),  # brewer.pal() 用 Set3 调色板为各通路分配连带颜色；ncol() 返回列数。
        process.label = 8,  # 将通路图例标签大小设为 8。
        border.size = 0)  # 去掉连带边框。

# 九象限图----------------------------------------------------------------------
#library("BiocManager")  # 备用代码，未执行；library() 加载已安装的包。
#BiocManager::install("dplyr")  # 备用代码，未执行；BiocManager::install() 安装或更新软件包。
#BiocManager::install("ggrepel")  # 备用代码，未执行；BiocManager::install() 安装或更新软件包。

# ---- 59. 两组学输入与对齐 ----
# 目的：按共同标识合并 RNA 和蛋白效应表，比较方向须一致。
library(dplyr)  # library() 加载 dplyr，用于整理、筛选和分组数据。
library(ggplot2)  # library() 加载 ggplot2，用于按图层绘图。
library(ggrepel)  # library() 加载 ggrepel，用于添加自动避让的文字标签。
conflicts_prefer(clusterProfiler::filter)  # conflicts_prefer() 优先选 clusterProfiler::filter，需确认其筛选行为；clusterProfiler::filter() 按逻辑条件筛选行。

pro <- read.table("8.九象限图/Protein.txt",header = T,sep = "\t")  # read.table() 读取蛋白 log2FC 与 P 值表，使用制表符分隔。
rna <- read.table("8.九象限图/RNA-seq.txt",header = T,sep = "\t")  # read.table() 读取 RNA log2FC 与 P 值表，比较方向须与蛋白一致。
head(pro,3)  # head() 查看数据前3行，核对结构和内容。
head(rna,3)  # head() 查看数据前3行，核对结构和内容。
combine <- merge(rna,pro,by.x="row",by.y="row",suffixes = c("_RNA","_Protein"))  # merge() 按 row 内连接两组学，给同名字段加 RNA/Protein 后缀；c() 将多个值合成向量。
dim(combine)  # dim() 查看数据行数和列数，核对输入规模。

# 区分象限

# ---- 60. 九象限幅度分类 ----
# 目的：以两组学 log2FC=±1 为界，将九个位置合并为五类变化模式。
combine <- combine %>%  # 将合并数据传入下一步，准备划分象限。
  mutate(  # mutate() 开始添加象限分类列 part。
    part = case_when(fc_RNA < -1 & fc_Protein > 1 ~ "part1",  # case_when() RNA 下调、蛋白上调，归为反向变化 part1。
                     between(fc_RNA, -1, 1) & fc_Protein > 1 ~ "part2",  # between() RNA 在 ±1 内、蛋白大于 1，归为 part2。
                     fc_RNA > 1 & fc_Protein > 1 ~ "part3",  # RNA 和蛋白均大于 1，归为同向变化 part3。
                     fc_RNA < -1 & between(fc_Protein, -1, 1) ~ "part4",  # between() 判断是否处于含端点的区间；RNA 小于 −1、蛋白在 ±1 内，归为 part4。
                     between(fc_RNA, -1, 1) & between(fc_Protein, -1, 1) ~ "part5",  # between() 两组学 log2FC 均在 ±1 内，归为小幅变化 part5。
                     fc_RNA > 1 & between(fc_Protein, -1, 1) ~ "part4",  # between() 判断是否处于含端点的区间；RNA 大于 1、蛋白在 ±1 内，同样归为 part4。
                     fc_RNA < -1 & fc_Protein < -1 ~ "part3",  # RNA 和蛋白均小于 −1，归为同向变化 part3。
                     between(fc_RNA, -1, 1) & fc_Protein < -1 ~ "part2",  # between() RNA 在 ±1 内、蛋白小于 −1，同样归为 part2。
                     fc_RNA > 1 & fc_Protein < -1 ~ "part1"  # RNA 上调、蛋白下调，归为反向变化 part1。
                     #TRUE ~ "#FFFFFF" # 默认白色
    )  # 结束条件分类；未匹配的记录会得到 NA。
  )  # 完成 part 列并保存到 combine。


# ---- 61. 两组学 P 值筛选 ----
# 目的：区分两个原始 P 值均达阈值的成员与其余成员。
table(combine$part)  # table() 统计各象限类别的成员数量。
head(combine,3)  # head() 查看数据前3行，核对结构和内容。
data <- combine[,c(1:3,5,6,8)]  # c() 将多个值合成向量；按列号提取标识、两组学效应和 P 值，须核对列名。
head(data,3)  # head() 查看数据前3行，核对结构和内容。
data$sig <- case_when(data$pvalue_RNA < 0.05 & data$pvalue_Protein < 0.05 ~ "sig",  # case_when() 两组学原始 P 值均小于 0.05，标为 sig。
                      data$pvalue_RNA >= 0.05 | data$pvalue_Protein >=0.05 ~ "no")  # 任一原始 P 值不小于 0.05，标为 no。
table(data$sig)  # table() 统计 sig/no 两类数量。
sig <- filter(data,sig == "sig")  # filter() 提取两组学原始 P 值均达阈值的记录。
non <- filter(data,sig == "no")  # filter() 提取未同时通过两个 P 值阈值的记录。

# ---- 62. 九象限基础散点 ----
# 目的：在 RNA 与蛋白效应坐标上显示全部成员，并按分类着色。
p <-ggplot(data,aes(fc_RNA,fc_Protein))+  # ggplot() 横轴为 RNA log2FC，纵轴为蛋白 log2FC；aes() 将数据变量映射到坐标或外观。
  geom_point(data=non,aes(fc_RNA,fc_Protein,color = "Filtered"))+  # geom_point() 将未同时通过 P 值阈值的点标为 Filtered；aes() 将数据变量映射到坐标或外观。
  geom_point(data=sig,aes(fc_RNA,fc_Protein,color=part),size=1.5)  # geom_point() 对双组学均达 P 值阈值的点按 part 着色；aes() 将数据变量映射到坐标或外观。
#guides(color="none") #图例
p  # 显示图形对象 p；在脚本中运行时可显式 print()。

# ---- 63. 九象限主题与边界 ----
# 目的：统一类别颜色、边距和坐标范围，并添加 ±1 分界线。
mycolor <- c("gray90","#3994F5","#FAAB31","#E94549","#E47EE1","#808080")  # c() 将多个值合成向量；定义 Filtered 和五类 part 对应的颜色。
top.mar=0.2  # 将上边距设为 0.2 英寸。
right.mar=0.2  # 将右边距设为 0.2 英寸。
bottom.mar=0.2  # 将下边距设为 0.2 英寸。
left.mar=0.2  # 将左边距设为 0.2 英寸。
mytheme<-theme_bw()+theme(text=element_text(family = "sans",colour ="gray30",size = 12),  # theme_bw() 使用白底主题，基础字体为 sans、深灰色、12 号；theme() 调整图形外观；element_text() 设置主题文字样式。
                          panel.grid = element_blank(),  # element_blank() 去掉绘图区的全部网格线。
                          panel.border = element_rect(size = 0.8,colour = "gray30"),  # element_rect() 添加宽度为 0.8 的深灰色面板边框。
                          axis.line = element_blank(),axis.ticks = element_line(size = 0.6,colour = "gray30"),  # element_blank() 隐藏轴基线，将刻度线设为宽 0.6 的深灰色；element_line() 设置主题线条样式。
                          axis.ticks.length = unit(1.5,units ="mm"),  # unit() 给长度附加单位；坐标轴刻度短线长度设为 1.5 mm。
                          plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),units="inches"))  # unit() 给长度附加单位；c() 将多个值合成向量；按上、右、下、左顺序应用英寸单位的边距。
p1 <- p + scale_colour_manual(name="",values=alpha(mycolor,0.7))+mytheme  # scale_colour_manual() 手动指定分类颜色；alpha() 给颜色设置透明度；为原图应用自定义主题和透明度为 0.7 的分类配色。
p1  # 显示图形对象 p1；在脚本中运行时可显式 print()。
p2 <- p1+geom_hline(yintercept = c(-1,1),size=0.5,color="gray40",lty="dashed")+  # geom_hline() 添加水平参考线；c() 将多个值合成向量；在蛋白 log2FC=±1 处添加水平虚线。
  geom_vline(xintercept = c(-1,1),size=0.5,color="gray40",lty="dashed")+  # geom_vline() 在 RNA log2FC=±1 处添加竖直虚线；c() 将多个值合成向量。
  scale_y_continuous(expand = expansion(add = c(0.5,0.5)),limits = c(-15,15),breaks = c(-15,-10,-5,0,5,10,15),label = c("-15","-10","-5","0","5","10","15"))+  # scale_y_continuous() 纵轴范围设为 −15 至 15，按指定数值标注刻度；expansion() 设置坐标两端留白；c() 将多个值合成向量。
  scale_x_continuous(expand = expansion(add = c(0.5,0.5)),limits = c(-15,15),breaks = c(-15,-10,-5,0,5,10,15),label = c("-15","-10","-5","0","5","10","15"))+  # scale_x_continuous() 横轴范围设为 −15 至 15，按指定数值标注刻度；expansion() 设置坐标两端留白；c() 将多个值合成向量。

# ---- 64. 标注变化较大成员 ----
# 目的：标出两组学 P 值均较小且至少一组学变化较大的成员。
  geom_text_repel(  # geom_text_repel() 添加自动避让的基因标签图层。
    data = subset(data, data$pvalue_RNA < 0.001 & data$pvalue_Protein < 0.001 & (abs(fc_RNA) > 5 | abs(fc_Protein) > 5)),  # subset() 标注两者 P<0.001 且至少一者 |log2FC|>5 的成员；abs() 取绝对值。
    aes(label = row),  # aes() 标签内容取 row 列中的成员标识。
    size = 2,  # 标签字号设为 2。
    color = "black",  # 将标签文字设为黑色。
    box.padding = 0.5,  # 标签框四周留白设为 0.5。
    max.overlaps = 50,  # 标签允许的重叠数量上限设为 50。
    segment.color = "grey60",  # 将标签连线设为灰色。
    min.segment.length = 0.2  # 仅绘制长度达到 0.2 的标签连线。
  )  # 结束本次函数调用。
p2  # 显示图形对象 p2；在脚本中运行时可显式 print()。

# ---- 65. 效应相关与图例 ----
# 目的：计算全部匹配成员的效应相关，并设置五类变化模式的显示标签。
cor = round(cor(data$fc_RNA,data$fc_Protein),2)  # round() 计算全部匹配成员的 Pearson 相关，保留两位小数；cor() 计算列向量间的相关系数。
lab <- paste("correlation=",cor,sep = "")  # paste() 将相关系数拼接成图中显示的文字。
lab  # 查看相关系数标签内容。
p3 <- p2+geom_text(x = -10, y = 10, label = lab, color="gray40")+  # geom_text() 添加文字图层；在坐标 (−10,10) 添加相关系数标签。
  scale_color_manual(  # scale_color_manual() 手动设置分类颜色和图例标签。
    name = "",  # 隐藏颜色图例标题。
    values = c("gray90","#3994F5","#FAAB31","#E94549","#E47EE1","#808080"),  # c() 将多个值合成向量；按 Filtered、part1 至 part5 的顺序分配颜色。
    labels = c("Filtered",  # c() 将多个值合成向量；第一类图例显示为 Filtered。
               "Reversed",  # part1 图例显示为 Reversed，表示反向变化。
               "Single_Sig",  # part2 显示为 Single_Sig，实际仅蛋白幅度越过 ±1。
               "Same trend",  # part3 显示为 Same trend，表示同向变化。
               "Single_Sig",  # part4 显示为 Single_Sig，实际仅 RNA 幅度越过 ±1。
               "none"  # part5 显示为 none，表示两者幅度均未越过 ±1。
    )  # 结束本次函数调用。
  ) +  # 结束本次设置，继续添加图层。
  theme(  # theme() 继续设置图形主题。
    legend.position = c(0.88, 0.18),  # c() 将多个值合成向量；将图例定位在图内相对坐标 (0.88,0.18)。
    legend.background = element_rect(fill = NA),  # element_rect() 将图例整体背景设为透明。
    legend.key = element_rect(fill = NA),  # element_rect() 将图例色块背景设为透明。
    legend.text = element_text(size = 8, color = "black",family = "sans")  # element_text() 图例标签设为黑色 sans 字体，字号为 8。
  )  # 结束本次函数调用。


# ---- 66. 保存九象限图 ----
# 目的：输出包含象限边界、成员标签和效应相关的图形。
pdf("九象限图_Gene.pdf",width = 6,height = 6)  # pdf() 打开 6×6 英寸的 PDF 图形设备。
png("九象限图_Gene.png",width = 6,height = 6,units = "in",res = 300)  # png() 打开 6×6 英寸、300 dpi 的 PNG 图形设备。
p3  # 显示图形对象 p3；在脚本中运行时可显式 print()。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

# PPI网络分析####

# ---- 67. PPI 数据和图对象 ----
# 目的：读取 STRING 边表，转换为图对象并计算节点连接度。
library(tidyverse)  # library() 加载 tidyverse，用于整理数据和绘图。
library(Biostrings)  # library() 加载 Biostrings，用于处理生物序列，此处为预加载。
library(ggfun)  # library() 加载 ggfun，用于提供绘图辅助工具。
library(ggraph)  # library() 加载 ggraph，用于绘制网络节点与连接。
library(tidygraph)  # library() 加载 tidygraph，用于整理图对象并计算节点属性。
library(purrr)  # library() 加载 purrr，用于批量处理列表。
library(rlang)  # library() 加载 rlang，用于提供表达式和函数编程工具。
library(patchwork)  # library() 加载 patchwork，用于组合多个图形。

# 读取从 STRING 得到的蛋白互作数据
links <- read.delim("9.ppi网络分析/string_interactions_short.tsv")  # read.delim() 读取 STRING 导出的蛋白相互作用边表。

graph_df <- as_tbl_graph(links) %>%  # as_tbl_graph() 将边表转为图对象，传入下一步计算节点属性。
  tidygraph::mutate(Popularity = centrality_degree(mode = 'out'))  # tidygraph::mutate() 计算节点出度并存为 Popularity；出度不代表生物调控方向；centrality_degree() 计算节点连接度。

# 保存网络图为 PDF 文件

# ---- 68. PPI 环形布局与节点 ----
# 目的：用圆形布局展示互作网络，以连接度表示节点大小。
pdf(file="PPI_1.pdf",width = 5.5,height = 5)  # pdf() 打开 5.5×5 英寸的 PDF 图形设备。
png("PPI_1.png",width = 5.5,height = 5,units = "in",res = 300)  # png() 打开 5.5×5 英寸、300 dpi 的 PNG 图形设备。
p1 <- ggraph(graph_df, layout = 'linear', circular = TRUE) +  # ggraph() 使用线性排列后围成圆形的网络布局。
  geom_edge_arc(aes(color = combined_score), alpha = 0.5) +  # geom_edge_arc() 绘制弧形边，按 STRING 综合置信分数着色，透明度为 0.5；aes() 将数据变量映射到坐标或外观。
  geom_node_point(aes(size = Popularity), fill= "#9e9ac8", color = "#000000", alpha=1, shape = 21) +  # geom_node_point() 节点大小映射出度，填充紫色并添加黑色边框；aes() 将数据变量映射到坐标或外观。
  scale_edge_color_continuous(low="#ffeda0", high="#e31a1c")+  # scale_edge_color_continuous() 边颜色从低分浅黄渐变到高分红色。
  #scale_edge_colour_distiller(palette = "RdPu", direction = 1) +   # 备用代码，未执行；scale_edge_colour_distiller() 使用连续调色板设置网络边颜色。
  coord_fixed() +  # coord_fixed() 固定横纵轴比例，避免圆形网络拉伸变形。
  scale_size(range = c(2.5,8)) +  # scale_size() 将节点大小的显示范围设为 2.5 至 8；c() 将多个值合成向量。

# ---- 69. PPI 标签与画布 ----
# 目的：让节点标签向外排列，并调整画布以容纳名称。
  geom_node_text(  # geom_node_text() 添加网络节点文字标签。
    aes(  # aes() 开始设置节点标签的内容、角度和位置映射。
      label = name,  # 标签内容取节点 name。
      angle = (atan2(y, x) * 180/pi) %% 180,  # atan2() 根据横纵坐标计算方位角；根据节点坐标计算标签角度，并限制在 0 至 180 度。
      x = x * 1.1,  # 将标签横坐标向外放大为节点坐标的 1.1 倍。
      y = y * 1.1,  # 将标签纵坐标向外放大为节点坐标的 1.1 倍。
    ),  # 结束本次函数调用。
    size = 3,  # 节点标签字号设为 3。
    hjust = "outward",  # 横向向外对齐标签。
    vjust = "outward",  # 纵向向外对齐标签。
    family = "serif",  # 使用 serif 衬线字体。
    check_overlap = FALSE  # 不自动隐藏发生重叠的标签。
  ) +  # 结束本次设置，继续添加图层。
  theme_nothing() +  # theme_nothing() 使用空白主题，去掉坐标轴等背景元素。
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))  # expand_limits() 将横纵坐标显示范围至少扩展到 ±1.3；c() 将多个值合成向量。
p1  # 显示图形对象 p1；在脚本中运行时可显式 print()。
dev.off()  # dev.off() 关闭当前图形设备，完成文件写出。

# kegg通路可视化
#BiocManager::install("pathview")  # 备用代码，未执行；BiocManager::install() 安装或更新软件包。

# ---- 70. 通路着色输入 ----
# 目的：读取蛋白表达并转换为猪 ENTREZID，连接注释与样本数值。
library(pathview)  # library() 加载 pathview，用于将数值映射到 KEGG 通路图。
library(tibble)  # library() 加载 tibble，用于转换行名与普通列。
library(clusterProfiler)  # library() 加载 clusterProfiler，用于进行标识转换和功能富集。
library(org.Ss.eg.db)  # library() 加载 org.Ss.eg.db，用于查询猪基因注释。

#注意读入的数据需要是归一化后的表达量
pathview_pro <- read.csv("10.KEGG通路可视化/kegg_pro.csv",header = T)  # read.csv() 读取用于 KEGG 通路着色的蛋白表达表。
pathview_pro$pro.Name <- toupper(pathview_pro$pro.Name)  # toupper() 将蛋白对应基因符号统一转为大写。
Pro <- bitr(pathview_pro$pro.Name,  # bitr() 开始将基因符号映射为猪 ENTREZID。
            fromType = "SYMBOL",  # 输入标识为基因符号 SYMBOL。
            toType = c("ENTREZID"),  # c() 将多个值合成向量；将基因符号转换为 ENTREZID。
            OrgDb = org.Ss.eg.db)  # 使用猪注释数据库 org.Ss.eg.db。
colnames(Pro) <- c("pro.Name","ENTREZID")  # colnames() 读取或设置列名；c() 将多个值合成向量；将映射表两列命名为 pro.Name 和 ENTREZID。
merged_df <- merge(Pro, pathview_pro, by = "pro.Name", all = TRUE)  # merge() 按 pro.Name 全连接映射与表达表，未匹配部分保留 NA。

library(conflicted)  # library() 加载 conflicted，用于指定同名函数优先使用的来源。

# ---- 71. 清理通路输入 ----
# 目的：去除未映射和重复 ID，仅保留带明确成员身份的表达矩阵。
conflicts_prefer(tidygraph::filter)  # conflicts_prefer() 优先选 tidygraph::filter，需确认其适用于当前数据框；tidygraph::filter() 按逻辑条件筛选行。
clean_merge <- merged_df %>%  # 将合并表传入清理步骤。
  filter(!is.na(ENTREZID)) %>%  # filter() 删除 ENTREZID 缺失的记录；is.na() 检查缺失值。
  mutate(ENTREZID = as.character(ENTREZID))  # mutate() 将 ENTREZID 转为字符型，便于按标识匹配；as.character() 转为字符型。

clean_merge <- clean_merge[!duplicated(clean_merge$ENTREZID),]  # duplicated() 标记重复记录；按 ENTREZID 去重，保留第一条记录。
rownames(clean_merge) = clean_merge$ENTREZID  # rownames() 读取或设置行名；将 ENTREZID 设为表达表行名。
clean_merge1 = clean_merge[,-c(1:2)]  # c() 将多个值合成向量；删除前两列标识信息，仅保留样本表达列。


# ---- 72. 数值标准化 ----
# 目的：按样本列跨蛋白做 z-score，明确后续颜色表示标准化样本值。
scaled <- scale(clean_merge1)  # scale() 默认按列减均值并除以标准差；对每个样本列跨蛋白做 z-score；此处不是计算组间 logFC。


# ---- 73. 多样本 KEGG 通路图 ----
# 目的：将各样本标准化值映射到猪糖酵解/糖异生通路节点。
pathview(gene.data = scaled,  # pathview() 将标准化样本值映射到 KEGG 通路图节点。
         pathway.id = "00010",  # 选择通路 00010，即糖酵解/糖异生通路。
         species = "ssc",  # 指定物种为猪 ssc。
         limit = list(gene = c(-1, 1)),  # list() 将对象组成列表；c() 将多个值合成向量；将基因数值着色范围设为 −1 至 1。
         low = "green",  # 将低值映射为绿色。
         mid = "gray",  # 将中间值映射为灰色。
         high = "red",  # 将高值映射为红色。
         node.sum = "mean")  # 同一节点对应多个输入成员时，取均值汇总。


# ---- 74. 单样本原生通路图 ----
# 目的：用第一个样本的数值绘制并保存 KEGG 原生图。
p <- pathview(gene.data = scaled[, 1], pathway.id = "00010", species = "ssc",  # pathview() 用第一个样本的标准化表达值绘制猪 00010 通路。
              out.suffix = "gse16873", kegg.native = T)  # 使用 KEGG 原生图，输出文件后缀设为 gse16873。

# 以分离图层方式再输出一个通路图版本；不是自动换成蛋白符号。

# ---- 75. 单样本分层通路图 ----
# 目的：将注释绘制为独立图层，另存同一通路的第二个版本。
p <- pathview(gene.data = scaled[, 1], pathway.id = "00010", species = "ssc",  # pathview() 再用第一个样本的数值绘制同一猪通路。
              out.suffix = "gse16873.2layer", kegg.native = T, same.layer = F)  # 设置新后缀并分层绘制；same.layer=FALSE 不负责转换名称。

