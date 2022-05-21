rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
setwd("D:/R/R.Projects/TCGA")
library(GEOquery)
gset <- getGEO('GSE154425', destdir=".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F) 
count1 <- exprs(gset[[1]])     ###表达矩阵文件

row.names(count1) <- count1[,1]   ##取行名
count1 <- count1[,-1]
count1 <- log2(count1+1)
##分组

pd=pData(gset[[1]])
group_list <- (pd$`treatment:ch1`)
                group_list <- as.factor(substr(colnames(tpms),1,3))
###   取编码蛋白的基因进行差异分析
Ginfo_0 <- read.table(file = "d:/R/R.Projects/gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo_0 <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] 
table(duplicated(Ginfo_0$genename))
Ginfo_0 <- Ginfo_0[!duplicated(Ginfo_0$genename),]
row.names(Ginfo_0) <- Ginfo_0$genename
comgene <- intersect(rownames(count1),rownames(Ginfo_0))
count1 <- count1[comgene,]


##根据列数判断TURE（=1）的个数，去掉大部分值为0的行
a[(apply(a,MARGIN = 1,function(x) sum(x>1)>3)),]   
table(apply(a,MARGIN = 1,function(x) sum(x>1)>3))
  a <- a[(apply(a,MARGIN = 1,function(x) sum(x>1)>3)),]  
library(org.Hs.eg.db)
ensembl <- toTable(org.Hs.egENSEMBL)
id <- toTable(org.Hs.egSYMBOL)
colnames(ensembl)
colnames(id)
##下面需要查看a与ensembl是否有相同列及列名
count1[,7] <- rownames(count1)
colnames(count1)[7] <- "ensembl_id"
merge(count1,ensembl,"ensembl_id",all.x=T)
tem <- merge(count1,ensembl,"ensembl_id",all.x=T)  ##需要相同的列名才能按照相同行把列合并
tem <- merge(tem,id,"gene_id",all.x=T)
##去除重复
table(duplicated(tem$symbol))
tem <- tem[!duplicated(tem$symbol),]    ##取剩下不重复的数据框
tem <- tem[!is.na(tem$symbol),]

##变为只有数值的数据框
row.names(tem) <- tem[,8]
tem <- tem[,-c(1,8)]

## hclust 
colnames(exprSet)=paste(group_list,1:18,sep='')
# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
hc=hclust(dist(t(exprSet)))
par(mar=c(5,5,5,10)) 
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)

## PCA 检测数据来源及分组可靠性
library(ggplot2)
library(ggfortify)
exprSet <- tem
df=as.data.frame(t(exprSet))
df$group=group_list 
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 
df=as.data.frame(t(exprSet))
dat.pca <- PCA(df, graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)


# DEG by limma 
suppressMessages(library(limma)) 
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(tem)
design
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix<-makeContrasts("era-ctr",levels = design)
contrast.matrix 
##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
##step1
fit <- lmFit(tem,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)


save.image(file = "note.Rdata")

