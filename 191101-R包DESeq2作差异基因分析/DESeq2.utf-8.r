#安装 DESeq2，常规方法 bioconductor 安装
BiocManager::install('DESeq2')

#数据量大时，嫌运行速率慢？试试最新版（v1.25.9），效率提升了数十倍
#但是最新版尚未添加至 bioconductor 里（bioconductor 安装的还是旧版），需要在 GitHub 中获取，如下安装
devtools::install_github('mikelove/DESeq2@ae7c6bd')

#若中间提示有其它依赖 R 包的旧版包冲突的话，先删除旧包再安装新的
remova.packages('xxx')
BiocManager::install('xxx')

###############################################
library(DESeq2)

#基因表达矩阵
gene <- read.delim('gene.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#指定分组因子顺序
#注意要保证表达矩阵中的样本顺序和这里的分组顺序是一一对应的
coldata <- data.frame(group = factor(rep(c('control', 'treat'), each = 8), levels = c('control', 'treat')))

##DESeq2 默认流程
#第一步，构建 DESeqDataSet 对象，详见 ?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = gene, colData = coldata, design = ~group)

#查看归一化后的 count 值分布
boxplot(log10(assays(dds)[['cooks']]), range = 0, las = 2)
plotDispEsts(dds)

#获取归一化的基因表达矩阵
vsd <- assay(vst(dds, blind = FALSE))
head(vsd, 10)
#write.table(vsd, 'norm_matrix.txt', sep = '\t', col.names = NA, quote = FALSE)

#第二步，差异分析，详见 ?DESeq 和 ?results
#标准方法
dds <- DESeq(dds, parallel = FALSE)	#parallel = TRUE 将启用多线程模式
suppressMessages(dds)

res <- results(dds, contrast = c('group', 'treat', 'control'), pAdjustMethod = 'fdr', alpha = 0.05)

#an alternate analysis: likelihood ratio test
ddsLRT <- DESeq(dds, test = 'LRT', reduced = ~ 1)
suppressMessages(ddsLRT)

resLRT <- results(ddsLRT, contrast = c('group', 'treat', 'control'), pAdjustMethod = 'fdr', alpha = 0.05)

#简要查看结果，例如
res
#再如
summary(res)
#再如
plotMA(res, alpha = 0.05, ylim = c(-3, 3))

#可以先按校正和 p 值由小到大排个序，方便查看
deseq_res <- as.data.frame(res[order(res$padj), ])

#输出
deseq_res$gene_id <- rownames(deseq_res)
write.table(deseq_res[c(7, 1:6)], 'DESeq2.txt', row.names = FALSE, sep = '\t', quote = FALSE)

##ggplot2 差异火山图
library(ggplot2)

deseq_res <- read.delim('DESeq2.txt', sep = '\t')

#例如这里根据 |log2FC| >= 1 & FDR p-value < 0.05 定义“差异”
deseq_res[which(deseq_res$padj %in% NA),'sig'] <- 'no diff'
deseq_res[which(deseq_res$log2FoldChange >= 1 & deseq_res$padj < 0.05),'sig'] <- 'rich (p.adj < 0.05, log2FC >= 1)'
deseq_res[which(deseq_res$log2FoldChange <= -1 & deseq_res$padj < 0.05),'sig'] <- 'down (p.adj < 0.05, log2FC <= -1)'
deseq_res[which(abs(deseq_res$log2FoldChange) < 1 | deseq_res$padj >= 0.05),'sig'] <- 'no diff'

#纵轴为显著性 p 值
volcano_p <- ggplot(deseq_res, aes(log2FoldChange, -log(padj, 10))) +
geom_point(aes(color = sig), alpha = 0.6, size = 1) +
scale_color_manual(values = c('blue2', 'gray30', 'red2')) +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.position = c(0.26, 0.92)) +
theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) + 
geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
labs(x = 'log2 Fold Change', y = '-log10 p-value', color = NA) +
xlim(-5, 5)

#ggsave('volcano_p.pdf', volcano_p, width = 5, height = 6)
ggsave('volcano_p.png', volcano_p, width = 5, height = 6)

#纵轴为基因表达值的 log10
volcano_count <- ggplot(deseq_res, aes(log2FoldChange, log10(baseMean))) +
geom_point(aes(color = sig), alpha = 0.6, size = 1) +
scale_color_manual(values = c('blue2', 'gray30', 'red2')) +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.position = c(0.2, 0.9)) +
theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) + 
labs(x = 'log2 Fold Change', y = 'Average log10 baseMean') +
xlim(-5, 5) +
coord_flip()

#ggsave('volcano_count.df', volcano_count, width = 7, height = 5)
ggsave('volcano_count.png', volcano_count, width = 7, height = 5)
