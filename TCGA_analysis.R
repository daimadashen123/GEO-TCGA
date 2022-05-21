rm(list = ls()) 
setwd("D:/R/R.Projects/TCGA")
library(tidyverse)
count1 <- read.table(file = "d:/R/R.Projects/TCGA.BLCA.sampleMap_HiSeqV2.gz",sep = "",header = T,fill = T)
row.names(count1) <- count1[,1]
count1 <- count1[,-1]
count1 <- count1[(apply(count1,MARGIN = 1,function(x) sum(x>0)>204)),]    ###去除检测表达较低的基因
group_list <- ifelse(as.numeric(str_sub(colnames(count1),14,15))<10,"tumor","normal") 
group_list <- factor(group_list,levels = c("normal","tumor")) 
table(group_list)

colnames(count1) <- c(paste0("ctr",1:3),paste0("era",1:3))
group_list <- as.factor(substr(colnames(count1),1,3))
###   取编码蛋白的基因进行差异分析
Ginfo_0 <- read.table(file = "d:/R/R.Projects/gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo_0 <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] 
table(duplicated(Ginfo_0$genename))
Ginfo_0 <- Ginfo_0[!duplicated(Ginfo_0$genename),]
row.names(Ginfo_0) <- Ginfo_0$genename
comgene <- intersect(rownames(count1),rownames(Ginfo_0))
count1 <- count1[comgene,]

###    DESeq2以及EdgeR包、limma包要求的输入文件为【未经标准化的原始read counts表达矩阵】！
###    DESeq2会对原始counts进行标准化（该方法优于CPM/TPM/FPKM/RPKM等方法），
###    去除不同测序批次导致的批次效应（batch effect），因此不能使用FPKM等已经标准化过的矩阵。
###    如果拿不到原始counts则转换为TPM后用limma包做差异分析。。。以下是fpkm转换为TPM做limma包
expMatrix <- count1
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)
exprSet <- tpms
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
exprSet <- log2(exprSet+1)
tem <- exprSet
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


###   DESeq2方法做差异分析
library(DESeq2) 
count1 <- ceiling(2^(count1)-1)
colData <- data.frame(row.names=colnames(count1), condition=group_list) #condition是你的分组信息 
dds <- DESeqDataSetFromMatrix(countData = count1, colData = colData, design = ~ condition) 
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast = c("condition",rev(levels(group_list)))) 
resOrdered <- res[order(res$pvalue),] 
DEG <- as.data.frame(resOrdered)
head(DEG)
DEG <- na.omit(DEG)
logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) ) 
logFC_cutoff 
DEG$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange)> 
                      logFC_cutoff,ifelse(DEG$log2FoldChange >
                      logFC_cutoff ,'UP','DOWN'),'NOT') ) 
head(DEG)
table(DEG$change)

###  edgeR方法做差异分析
library(edgeR) 
dge <- DGEList(counts=count1,group=group_list) 
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
DEG2=topTags(fit2, n=nrow(count1)) 
DEG2=as.data.frame(DEG2) 
logFC_cutoff2 <- with(DEG2,mean(abs(logFC)) + 2*sd(abs(logFC)) ) 
DEG2$change = as.factor(  ifelse(DEG2$PValue < 0.05 & abs(DEG2$logFC) >
                                  logFC_cutoff2,  ifelse(DEG2$logFC >logFC_cutoff2 ,'UP','DOWN'),'NOT')  ) 
head(DEG2)
DESeq2_DEG <- DEG
table(DEG2$change) 
edgeR_DEG <- DEG2

###   limma-voom方法做差异分析     
library(limma) 
design <- model.matrix(~0+group_list) 
colnames(design)=levels(group_list) 
rownames(design)=colnames(count1) 
dge <- DGEList(counts=count1) 
dge <- calcNormFactors(dge) 
logCPM <- cpm(dge, log=TRUE, prior.count=3) 
v <- voom(dge,design, normalize="quantile") 
fit <- lmFit(v, design) 
constrasts = paste(rev(levels(group_list)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit3=contrasts.fit(fit,cont.matrix) 
fit3=eBayes(fit3) 
DEG3 = topTable(fit3, coef=constrasts, n=Inf) 
DEG3 = na.omit(DEG3) 
logFC_cutoff3 <- with(DEG3,mean(abs(logFC)) + 2*sd(abs(logFC)) ) 
DEG3$change = as.factor(  ifelse(DEG3$P.Value < 0.05 & abs(DEG3$logFC) >
                                 logFC_cutoff3,  ifelse(DEG3$logFC > logFC_cutoff3 ,'UP','DOWN'),'NOT')  ) 
head(DEG3)
limma_voom_DEG <- DEG3
table(DEG3$change) 

###    保存三种方法得到的矩阵
save(DESeq2_DEG,edgeR_DEG,limma_voom_DEG,group_list,file = "DEG.Rdata")
load(file = "DEG.Rdata")
###    热图
cg1 = rownames(DESeq2_DEG)[DESeq2_DEG$change !="NOT"]
library(pheatmap)
library(RColorBrewer)
color<- colorRampPalette(c('#436eee','white','#EE0000'))(100) #下面画第一个DESeq2矩阵的热图 
mat=count1[cg1,] 
n=t(scale(t(mat))) 
n[n>1]=1 
n[n< -1]= -1 
ac=data.frame(group=group_list) 
rownames(ac)=colnames(mat) 
ht1 <- pheatmap(n,show_rownames = F,show_colnames = F, cluster_rows = F,cluster_cols = T, annotation_col = ac,color=color)

###   火山图
library(EnhancedVolcano) 
library(airway) 
EnhancedVolcano(DESeq2_DEG, lab = rownames(DESeq2_DEG), x = 'log2FoldChange', y = 'pvalue', xlim = c(-8, 8),ylim = c(0, -log10(10e-30)),pointSize = 0.8, labSize = 2.0,title = 'DESeq2_DEG', pCutoff = 10e-20, FCcutoff = 2.5,  col=c('black', 'blue', 'green', 'red1'), colAlpha =0.5, legendLabels=c('NS','log2FoldChange','P value', 'P value & log2FoldChange'), legendPosition = 'right', legendLabSize = 10, legendIconSize = 3.0, )
