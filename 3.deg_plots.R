rm(list = ls())
cancer_type = "Pancreas"
load("Pancreas_expr_ph.Rdata")
load("Pancreas_DEG.Rdata")
library(ggplot2)
library(tinyarray)
dat = Pancreas_expr
pca.plot = draw_pca(dat,group_list);pca.plot
save(pca.plot,file = paste0(cancer_type,"pcaplot.Rdata"))

cg1 = rownames(DESeq2_DEG)[DESeq2_DEG$change !="NOT"]
cg2 = rownames(edgeR_DEG)[edgeR_DEG$change !="NOT"]
cg3 = rownames(limma_voom_DEG)[limma_voom_DEG$change !="NOT"]

h1 = draw_heatmap(dat[cg1,],group_list,scale_before = T)
h2 = draw_heatmap(dat[cg2,],group_list,scale_before = T)
h3 = draw_heatmap(dat[cg3,],group_list,scale_before = T)

m2d = function(x){
  mean(abs(x))+2*sd(abs(x))
}

v1 = draw_volcano(DESeq2_DEG,pkg = 1,logFC_cutoff = m2d(DESeq2_DEG$log2FoldChange))
v2 = draw_volcano(edgeR_DEG,pkg = 2,logFC_cutoff = m2d(edgeR_DEG$logFC))
v3 = draw_volcano(limma_voom_DEG,pkg = 3,logFC_cutoff = m2d(limma_voom_DEG$logFC))

library(patchwork)
(h1 + h2 + h3) / (v1 + v2 + v3) +plot_layout(guides = 'collect')

ggsave(paste0(cancer_type,"heat_vo.png"),width = 15,height = 10)



UP=function(df){
  rownames(df)[df$change=="UP"]
}
DOWN=function(df){
  rownames(df)[df$change=="DOWN"]
}

up = intersect(intersect(UP(DESeq2_DEG),UP(edgeR_DEG)),UP(limma_voom_DEG))
down = intersect(intersect(DOWN(DESeq2_DEG),DOWN(edgeR_DEG)),DOWN(limma_voom_DEG))
dat = Pancreas_expr
hp = draw_heatmap(dat[c(up,down),],group_list,scale_before = T)

#上调、下调基因分别画维恩图

up.plot <- draw_venn(list(Deseq2 = UP(DESeq2_DEG),
                          edgeR = UP(edgeR_DEG),
                          limma = UP(limma_voom_DEG)),
                     "UPgene"
)
down.plot <- draw_venn(list(Deseq2 = DOWN(DESeq2_DEG),
                           edgeR = DOWN(edgeR_DEG),
                           limma = DOWN(limma_voom_DEG)),
                       "DOWNgene"
)
#维恩图拼图，终于搞定

library(patchwork)
#up.plot + down.plot
# 就爱玩拼图
pca.plot + hp+up.plot +down.plot
ggsave(paste0(cancer_type,"heat_ve_pca.png"),width = 15,height = 10)
