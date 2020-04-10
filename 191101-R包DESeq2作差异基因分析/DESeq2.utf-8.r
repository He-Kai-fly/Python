#��װ DESeq2�����淽�� bioconductor ��װ
BiocManager::install('DESeq2')

#��������ʱ�����������������������°棨v1.25.9����Ч����������ʮ��
#�������°���δ����� bioconductor �bioconductor ��װ�Ļ��Ǿɰ棩����Ҫ�� GitHub �л�ȡ�����°�װ
devtools::install_github('mikelove/DESeq2@ae7c6bd')

#���м���ʾ���������� R ���ľɰ����ͻ�Ļ�����ɾ���ɰ��ٰ�װ�µ�
remova.packages('xxx')
BiocManager::install('xxx')

###############################################
library(DESeq2)

#���������
gene <- read.delim('gene.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#ָ����������˳��
#ע��Ҫ��֤�������е�����˳�������ķ���˳����һһ��Ӧ��
coldata <- data.frame(group = factor(rep(c('control', 'treat'), each = 8), levels = c('control', 'treat')))

##DESeq2 Ĭ������
#��һ�������� DESeqDataSet ������� ?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = gene, colData = coldata, design = ~group)

#�鿴��һ����� count ֵ�ֲ�
boxplot(log10(assays(dds)[['cooks']]), range = 0, las = 2)
plotDispEsts(dds)

#��ȡ��һ���Ļ��������
vsd <- assay(vst(dds, blind = FALSE))
head(vsd, 10)
#write.table(vsd, 'norm_matrix.txt', sep = '\t', col.names = NA, quote = FALSE)

#�ڶ����������������� ?DESeq �� ?results
#��׼����
dds <- DESeq(dds, parallel = FALSE)	#parallel = TRUE �����ö��߳�ģʽ
suppressMessages(dds)

res <- results(dds, contrast = c('group', 'treat', 'control'), pAdjustMethod = 'fdr', alpha = 0.05)

#an alternate analysis: likelihood ratio test
ddsLRT <- DESeq(dds, test = 'LRT', reduced = ~ 1)
suppressMessages(ddsLRT)

resLRT <- results(ddsLRT, contrast = c('group', 'treat', 'control'), pAdjustMethod = 'fdr', alpha = 0.05)

#��Ҫ�鿴���������
res
#����
summary(res)
#����
plotMA(res, alpha = 0.05, ylim = c(-3, 3))

#�����Ȱ�У���� p ֵ��С�����Ÿ��򣬷���鿴
deseq_res <- as.data.frame(res[order(res$padj), ])

#���
deseq_res$gene_id <- rownames(deseq_res)
write.table(deseq_res[c(7, 1:6)], 'DESeq2.txt', row.names = FALSE, sep = '\t', quote = FALSE)

##ggplot2 �����ɽͼ
library(ggplot2)

deseq_res <- read.delim('DESeq2.txt', sep = '\t')

#����������� |log2FC| >= 1 & FDR p-value < 0.05 ���塰���족
deseq_res[which(deseq_res$padj %in% NA),'sig'] <- 'no diff'
deseq_res[which(deseq_res$log2FoldChange >= 1 & deseq_res$padj < 0.05),'sig'] <- 'rich (p.adj < 0.05, log2FC >= 1)'
deseq_res[which(deseq_res$log2FoldChange <= -1 & deseq_res$padj < 0.05),'sig'] <- 'down (p.adj < 0.05, log2FC <= -1)'
deseq_res[which(abs(deseq_res$log2FoldChange) < 1 | deseq_res$padj >= 0.05),'sig'] <- 'no diff'

#����Ϊ������ p ֵ
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

#����Ϊ������ֵ�� log10
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
