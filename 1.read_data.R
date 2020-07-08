rm(list = ls())
###1.读取全部数据----
if(F){
  counts = data.table::fread("TcgaTargetGtex_gene_expected_count.gz",data.table = F)
  ph = data.table::fread("TcgaTargetGTEX_phenotype.txt.gz",data.table = F)
  counts[1:4,1:4]
  rownames(counts) <- counts$sample
  counts = counts[,-1]
  rownames(ph) <- ph$sample
  ph = ph[,-1]
  expr = as.matrix(counts)
  dim(expr)
  dim(ph)
  #匹配样本名
  table(rownames(ph)%in%colnames(expr))
  table(colnames(expr)%in%rownames(ph))
  ph = ph[match(colnames(expr),rownames(ph)),]
  identical(colnames(expr),rownames(ph))
  #挑选出TCGA和gtex样本
  table(ph$`_study`)
  keep = ph$`_study`!="TARGET"
  ph=ph[keep,]
  expr = expr[,keep]
  save(expr,ph,file = "gtex_tcga_expr_ph.Rdata")
}
load("gtex_tcga_expr_ph.Rdata")
### 挑出胰腺癌数据
table(ph$`_primary_site`)
keep2 = ph$`_primary_site`=="Pancreas"
Pancreas_expr = expr[,keep2]
Pancreas_ph = ph[keep2,]
save(Pancreas_ph,Pancreas_expr,file = "Pancreas_expr_ph.Rdata")
