rm(list = ls()) 
setwd("D:/R/R.Projects/TCGA")

counttumor <- count1[,as.numeric(str_sub(colnames(count1),14,15))<10]
###********
library(GSVA)
library(GSEABase)
c2gmt <- getGmt("d:/R/R.Projects/c2.cp.v7.5.1.symbols.gmt")
gene.set <- c2gmt[grep("^KEGG|REACTOME|BIOCARTA", names(c2gmt)),]
gs.exp <- gsva(as.matrix(counttumor), gene.set, kcdf = "Poisson", min.sz = 10)
save(gs.exp,file = "gs.exp.Rdata")
####这里 kcdf 参数设为"Poisson"是因为使用read count数据，
####如果是使用 log 后的CPM, RPKM, TPM等数据就用默认值"Gaussian"。
cg1 = rownames(gs.exp)
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#436eee','white','#EE0000'))(100) #下面画第一个DESeq2矩阵的热图 
mat=count1[cg1,] 
n=t(scale(t(gs.exp))) 
n[n>2]=2 
n[n< -2]= -2 
ac=data.frame(group=group_list) 
rownames(ac)=colnames(mat) 
ht1 <- pheatmap(n,show_rownames = F,show_colnames = F, cluster_rows =T,cluster_cols = T,)








DEA.gs <- TCGAanalyze_DEA(
  mat1 = gs.exp[, colnames(luad.exp)],
  metadata = FALSE,
  pipeline = "limma",
  Cond1type = "LUAD",
  fdr.cut = 0.05,
  logFC.cut = 0.5,
)
#测试数据的生成参考了参考资料的第三个

