#差异分析----
###1.数据整理------------------
rm(list = ls())
load("Pancreas_expr_ph.Rdata")
#log逆转
Pancreas_expr = 2^Pancreas_expr -1
Pancreas_expr[1:4,1:4]
#tumor和normal统计
library(stringr)
k1 = str_starts(rownames(Pancreas_ph),"TCGA")
k2 = as.numeric(str_sub(rownames(Pancreas_ph),14,15))<10
table(k1&k2)
group_list = ifelse(k1&k2,"tumor","normal")
group_list = factor(group_list,levels = c("normal","tumor"))
Pancreas_expr[1:4,1:4]
k3 = apply(Pancreas_expr,1,function(x){sum(x>1)})>150
table(k3)
Pancreas_expr = Pancreas_expr[k3,]
###2.DESeq2---------
expr = floor(Pancreas_expr)
expr[1:4,1:4]
library(DESeq2)
colData <- data.frame(row.names =colnames(expr), 
                      condition=group_list)
dds <- DESeqDataSetFromMatrix(
  countData = expr,
  colData = colData,
  design = ~ condition)
#参考因子应该是对照组 dds$condition <- relevel(dds$condition, ref = "untrt")

dds <- DESeq(dds)
# 两两比较
res <- results(dds, contrast = c("condition",rev(levels(group_list))))
resOrdered <- res[order(res$pvalue),] # 按照P值排序
DEG <- as.data.frame(resOrdered)
head(DEG)
# 去除NA值
DEG <- na.omit(DEG)

#添加change列标记基因上调下调
logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
#logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
         ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)
head(DEG)
table(DEG$change)
DESeq2_DEG <- DEG
###3.edgeR---------
expr = Pancreas_expr
library(edgeR)

dge <- DGEList(counts=expr,group=group_list)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge) 

design <- model.matrix(~0+group_list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group_list)

dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast=c(-1,1)) 

DEG=topTags(fit2, n=nrow(exp))
DEG=as.data.frame(DEG)
logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
#logFC_cutoff <- 2
DEG$change = as.factor(
  ifelse(DEG$PValue < 0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
head(DEG)
table(DEG$change)
edgeR_DEG <- DEG
###limma----
library(limma)

design <- model.matrix(~0+group_list)
colnames(design)=levels(group_list)
rownames(design)=colnames(exp)

dge <- DGEList(counts=expr)
dge <- calcNormFactors(dge)

v <- voom(dge,design, normalize="quantile")
fit <- lmFit(v, design)

constrasts = paste(rev(levels(group_list)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)
logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
#logFC_cutoff <- 2
DEG$change = as.factor(
  ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
head(DEG)
limma_voom_DEG <- DEG
tj = data.frame(deseq2 = as.integer(table(DESeq2_DEG$change)),
                edgeR = as.integer(table(edgeR_DEG$change)),
                limma_voom = as.integer(table(limma_voom_DEG$change)),
                row.names = c("down","not","up")
);tj
save(DESeq2_DEG,edgeR_DEG,limma_voom_DEG,group_list,tj,file = "Pancreas_DEG.Rdata")
