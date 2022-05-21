
library(tidyverse)
https://xenabrowser.net/datapages/

##文件的读�?
#读取tsv文件
counts1 = read.table(file = 'TCGA-LIHC.htseq_counts.tsv', sep = '\t', header = TRUE) 
rownames(counts1) <- counts1[,1] #Alt <- 
x <- counts1[,1:3]
counts1 = counts1[,-1]
#substr函数
substr("wanglihong",1,4)
#table函数
table(substr(colnames(counts1),14,16))
#c("01A","11A")
#%in%符号用于判断是否属于
counts1 <- counts1[,substr(colnames(counts1),14,16)%in% c("01A","11A")]

table(substr(colnames(counts1),14,16))

#保留行名�?15�?
rownames(counts1) <- substr(rownames(counts1),1,15)
ceiling(1.2)
ceiling(3.8)
counts <- ceiling(2^(counts1)-1)

##文件的输�?
#输出为文�?
write.table(counts,"counts.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#输出为表�?
write.csv(counts, file = "counts.csv")

#9.11
####差异分析####
#设置工作目录
setwd("xena")
##读取文本文件
counts <- read.table("counts.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#加载基因注释文件
Ginfo_0 <- read.table(file = "d:/R/R.Projects/gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo_0 <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] #只要编码RNA
#美元符号代表提取�?
#取行名交�?
comgene <- intersect(rownames(counts),rownames(Ginfo))
counts <- counts[comgene,]
class(counts)#判断数据类型
class(comgene)
Ginfo <- Ginfo[comgene,]
a <- rownames(counts)
b <- rownames(Ginfo)
identical(a,b)

counts$Gene <- as.character(Ginfo$genename)   #新增Gene Symbol
counts <- counts[!duplicated(counts$Gene),]   #去重�?
rownames(counts) <- counts$Gene   #将行名变为Gene Symbol
ncol(Ginfo)
nrow
counts <- counts[,-ncol(counts)]   #去除最后一�?
write.table(counts, file = "LIHC_counts_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存癌症患者的counts
tumor <- colnames(counts)[substr(colnames(counts),14,16) == "01A"]
counts_01A <- counts[,tumor]
write.table(counts_01A, file = "LIHC_counts_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#差异分析
library(tidyverse)
#安装BiocManager
if(!require(DESeq2))BiocManager::install('DESeq2')
library(DESeq2)

counts = counts[apply(counts, 1, function(x) sum(x > 1) > 32), ]
conditions=data.frame(sample=colnames(counts),
                      group=factor(ifelse(substr(colnames(counts),14,16) == "01A","T","N"),levels = c("N","T"))) %>% 
  column_to_rownames("sample")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conditions,
  design = ~ group)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
save(res,file = "LIHC_DEG.rda")#一定要保存�?
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 3, padj < 0.05)#根据自己需�?
#DEG:differentially expressed genes
#读取差异基因文件
#直接在文件夹双击

####整理fpkm文件####
#与counts几乎相同，fpkm不需进行log转换
#读取tsv文件
library(tidyverse)
setwd("xena")
fpkm1 = read.table(file = 'TCGA-LIHC.htseq_fpkm.tsv', sep = '\t', header = TRUE) 
rownames(fpkm1) <- fpkm1[,1]  
fpkm1 = fpkm1[,-1]
table(substr(colnames(fpkm1),14,16))
fpkm1 <- fpkm1[,substr(colnames(fpkm1),14,16)%in% c("01A","11A")]
table(substr(colnames(fpkm1),14,16))
rownames(fpkm1) <- substr(rownames(fpkm1),1,15)
fpkm <- fpkm1

Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] #只要编码RNA
#取行名交�?
comgene <- intersect(rownames(fpkm),rownames(Ginfo))
fpkm <- fpkm[comgene,]
Ginfo <- Ginfo[comgene,]
fpkm$Gene <- as.character(Ginfo$genename)   #新增Gene Symbol
fpkm <- fpkm[!duplicated(fpkm$Gene),]   #去重�?
rownames(fpkm) <- fpkm$Gene   #将行名变为Gene Symbol
fpkm <- fpkm[,-ncol(fpkm)]   #去除最后一�?
#保存所以患者的fpkm文件
write.table(fpkm, file = "LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存癌症患者的fpkm文件
tumor <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "01A"]
fpkm_01A <- fpkm[,tumor]
write.table(fpkm_01A, file = "LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存正常样本的fpkm文件
normal <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "11A"]
fpkm_11A <- fpkm[,normal]
write.table(fpkm_11A, file = "LIHC_fpkm_mRNA_11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#整理完毕#

####9.14####
setwd("xena")
library(tidyverse)
fpkm_01A <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
fpkm_11A <- read.table("LIHC_fpkm_mRNA_11A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#读取之前的差异分析结果，还记得怎么读取�?#
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)

gene <- c("LIN28B","CTAG2","REG3A")
a <- fpkm_01A[gene,]
b <- fpkm_11A[gene,]
a <- t(a)
b <- t(b)
class(a)
a <- as.data.frame(a)
b <- as.data.frame(b)
##运用传导�?%>%  cltrl+shift+M 
a <- a %>% t() %>% as.data.frame()
b <- b %>% t() %>% as.data.frame()
write.csv(a, file = "01A.csv")
write.csv(b, file = "11A.csv")
#Graphpad等等


####GEO数据库的使用####
####代表什�?
#pubmed插件 文献检�?
#GEO网站：https://www.ncbi.nlm.nih.gov/geo/
####GSE84402####
setwd("GSE84402")
###加载R�?
library(tidyverse)
chooseBioCmirror()
BiocManager::install('GEOquery')
library(GEOquery)
###下载数据，如果文件夹中有会直接读�?
chooseBioCmirror()
gset = getGEO('GSE84402', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
###提取子集
gset[[1]]

#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)
library(stringr)
#设置参考水�?
group_list <- ifelse(str_detect(pdata$source_name_ch1, "hepatocellular carcinoma"), "tumor",
                     "normal")
#因子�?
group_list = factor(group_list,
                    levels = c("normal","tumor"))
##2.2 通过exprs函数获取表达矩阵并校�?
exp <- exprs(gset[[1]])
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
dev.off()
###数据校正
library(limma) 
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)
exp <- log2(exp+1)
range(exp)
dev.off()
#使用R包转换id
index = gset[[1]]@annotation
if(!require("hgu133plus2.db"))
  BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)
#length(unique(ids$symbol))
#table(sort(table(ids$symbol)))
#id转换
library(tidyverse)
exp <- as.data.frame(exp)
exp <- exp %>% mutate(probe_id=rownames(exp))
exp <- exp %>% inner_join(ids,by="probe_id") 
exp <- exp[!duplicated(exp$symbol),]
rownames(exp) <- exp$symbol
exp <- exp[,-(29:30)]
write.table(exp, file = "exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####GEO手动注释####
####GSE31056####
setwd("GSE31056")
###加载R�?
library(tidyverse)
BiocManager::install('GEOquery')
library(GEOquery)
###下载数据，如果文件夹中有会直接读�?
chooseBioCmirror()
gset = getGEO('GSE31056', destdir=".", AnnotGPL = F, getGPL = F)
#有时会报�?  Increase it by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
class(gset)
###提取子集
gset[[1]]
#读取表达�?
exp <- exprs(gset[[1]])
#把表达谱转为数据框格�?
exp <- as.data.frame(exp)
##转换id
#读取GPL文件
comname <- intersect(rownames(exp),rownames(GPL))
exp <- exp[comname,]
GPL <- GPL[comname,]
exp1 <- cbind(GPL,exp)
exp1 <- exp1[!duplicated(exp1$SYMBOL),]
rownames(exp1) <- exp1$SYMBOL
exp1 <- exp1[,-(1:5)]
write.table(exp1, file = "exp1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



####中秋快乐####
setwd("GSE84402")
###加载R�?
library(tidyverse)
library(GEOquery)
###下载数据，如果文件夹中有会直接读�?
gset = getGEO('GSE84402', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
###提取子集
gset[[1]]
#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)
library(stringr)
#设置参考水�?
group_list <- ifelse(str_detect(pdata$source_name_ch1, "hepatocellular carcinoma"), "tumor",
                     "normal")
#因子�?
group_list = factor(group_list,
                    levels = c("normal","tumor"))

##读取上节课整理好的表达数据exp##
exp <- read.table("exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#差异分析
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
write.table(deg, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
##标记上下调基�?
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)

##热图##
cg = rownames(deg)[deg$change !="stable"]
diff=exp[cg,]
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(diff) 
pheatmap(diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

####GEO三大富集分析
setwd("GSE84402")
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
deg <- read.table("deg_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
deg <- deg %>% filter(change!="stable")

DEG <- deg
DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

#GO分析
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result

#KEGG
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk_res <- kk@result

#GSEA
msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "c5.all.v7.0.entrez.gmt"    #c2.all.v7.0.entrez.gmt �? c5.all.v7.0.entrez.gmt
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

geneList = DEG[,2]
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)

set.seed(1)
KEGG<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
#转换成数据框
KEGG_result_df <- as.data.frame(KEGG)
write.table(KEGG_result_df,file="GSEA_MSigDb_C5_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
save(KEGG,KEGG_result_df,file = "GSEA_deg_SPP1.rda")


####cox回归分析####
#设置工作目录
setwd("cox")
#安装加载R�?
install.packages("survival")
install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)
#下载生存信息
#xena官网：https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Liver%20Cancer%20(LIHC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#读取生存信息tsv文件
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
#整理生存信息数据
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
#读取表达数据
expr <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]
#表达数据整理完毕
#读取tcga差异分析结果
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)
#整合
deg_expr <- expr[rownames(res_deseq2),] %>% t() %>% as.data.frame()
surv.expr <- cbind(surv,deg_expr)

#Cox分析
Coxoutput <- NULL 
for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}


write.table(Coxoutput, file = "cox results.txt",sep = "\t",row.names = F,col.names = T,quote = F)
###筛选top基因
pcutoff <- 0.001
topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),] # 取出p值小于阈值的基因
topgene <- topgene[1:10,]

#3. 绘制森林�?
##3.1 输入表格的制�?
tabletext <- cbind(c("Gene",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$p),3),nsmall = 3)))
##3.2 绘制森林�?
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(0.5,1,1.5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "12" = gpar(lwd=2, col="black")), # 在最后一行上画黑色实�?
           lineheight = unit(.75,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
dev.off()


####计算患者免疫评分与肿瘤纯度#####
setwd("TCGA ESTIMATE")  #设置工作目录
#安装�?
library(utils) #这个包应该不用下载，自带�?
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)
#读取肿瘤患�?01A表达�?
expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


#计算免疫评分
filterCommonGenes(input.f = "LIHC_fpkm_mRNA_01A.txt",   #输入文件�?
                  output.f = "LIHC_fpkm_mRNA_01A.gct",   #输出文件�?
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("LIHC_fpkm_mRNA_01A.gct",   #刚才的输出文件名
              "LIHC_fpkm_mRNA_01A_estimate_score.txt",   #新的输出文件名（即估计的结果文件�?
              platform="affymetrix")   #默认平台

#3. 输出每个样品的打�?
result <- read.table("LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
result <- result[,-1]   
colnames(result) <- result[1,]   
result <- as.data.frame(t(result[-1,]))

rownames(result) <- colnames(expr)
write.table(result, file = "LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得�?



####ROC####
#读取生存信息tsv文件
setwd("ROC")
library(tidyverse)
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
#整理生存信息数据
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
#保存整理好的生存信息
write.table(surv, file = "survival.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#读取表达数据
expr <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]
#提取上次cox作图�?10个基�?
gene <- c("SPP1","PAGE1","G6PD","MAGEA4",'CDCA8',
          'TRIM54','KIF2C','KIF20A','ANLN',"SLC7A11")
exp10 <- expr[gene,] %>% t() %>% as.data.frame()
#整合表达谱与生存信息
exp_sur <- cbind(exp10,surv)
write.table(exp_sur, file = "exp_sur.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#准备R�?
install.packages("ROCR")
install.packages("rms")
library(ROCR)
library(rms)
#构建ROC预测模型
ROC1 <- prediction(exp_sur$SPP1,exp_sur$OS)   #构建ROC预测模型 
ROC2 <- performance(ROC1,"tpr","fpr")   #计算预测模型的TPR/FPR�?
AUC <- performance(ROC1,"auc")   #计算曲线下面�?(AUC)�?

AUC<- 0.5604839 #�? 根据结果对AUC进行赋�?

#1.4 绘制ROC曲线
plot(ROC2,
     col="red",   #曲线的颜�?
     xlab="False positive rate", ylab="True positive rate",   #x轴和y轴的名称
     lty=1,lwd=3,
     main=paste("AUC=",AUC))
abline(0, 1, lty=2, lwd=3)   #绘制对角�?
dev.off()

####timeROC####
setwd("timeROC")
#R�?
install.packages("timeROC")
install.packages("survival")
library(timeROC)
library(survival)
library(tidyverse)

#2.2 数据的整理与载入
exp_sur <- read.table("exp_sur.txt", header=T,sep="\t", check.names=F, row.names=1)
exp_sur$OS.time <- exp_sur$OS.time/365
exp_sur_01A <- exp_sur[substr(rownames(exp_sur),14,16) == "01A",]
write.table(exp_sur_01A, file = "exp_sur_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#2.3 构建ROC曲线函数
ROC3 <- timeROC(T=exp_sur_01A$OS.time,   #结局时间
                delta=exp_sur_01A$OS,   #结局指标
                marker=exp_sur_01A$SPP1,   #预测变量
                cause=1,   #阳性结局指标数�?
                weighting="marginal",   #计算方法，默认为marginal
                times=c(1, 3, 5),   #时间点，选取1年，3年和5年的生存�?
                iid=TRUE)
ROC3   #查看模型变量信息

#2.4 绘制ROC曲线
plot(ROC3,
     time=1, col="red")   #time是时间点，col是线条颜�?
plot(ROC3,
     time=3, col="green", add=TRUE)   #add指是否添加在上一张图�?
plot(ROC3,
     time=5, col="blue", add=TRUE)
legend("bottomright",
       c("Year-1", "Year-3", "Year-5"),
       col=c("red", "green", "blue"),
       lty=1, lwd=2)   #添加标签信息

dev.off()



####TCGA差异分析热图####
setwd("xena")
library(tidyverse)
exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)

logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)
library(pheatmap)
cg = rownames(DEG)[DEG$change !="NOT"]
exp_diff <- exp[cg,]
group_list=factor(ifelse(substr(colnames(exp),14,16) == "01A","T","N"),levels = c("N","T"))
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(exp_diff)
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

####TCGA差异分析火山�?####
setwd("xena")
library(tidyverse)
exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)

logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)


install.packages("ggpubr")
install.packages("ggthemes")
library(ggpubr)
library(ggthemes)

DEG$logP <- -log10(DEG$padj)
ggscatter(DEG,
          x = "log2FoldChange", y = "logP") +
  theme_base()

#增加基因上下调信�?
ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base()

#添加分界�?
ggscatter(DEG, x = "log2FoldChange", y = "logP", xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")
dev.off()

#添加基因标签信息
DEG$Label = ""   #新加一列label
DEG <- DEG[order(DEG$padj), ]   #对差异基因的p值进行从小到大的排序
DEG$Gene <- rownames(DEG)
#高表达的基因中，选择fdr值最小的5�?
up.genes <- head(DEG$Gene[which(DEG$change == "UP")], 5)
#低表达的基因中，选择fdr值最小的5�?
down.genes <- head(DEG$Gene[which(DEG$change == "DOWN")], 5)
#将up.genes和down.genes合并，并加入到Label�?
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
DEG$Label[match(DEG.top5.genes, DEG$Gene)] <- DEG.top5.genes

ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1,
          label = DEG$Label,
          font.label = 8,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)") +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

dev.off()

####cibersort####
setwd("cibersort")   
install.packages('e1071')
install.packages('parallel')
#install.packages("BiocManager")
BiocManager::install("preprocessCore", version = "3.13")
library(e1071)
library(parallel)
library(preprocessCore)
source("CIBERSORT.R")   
sig_matrix <- "LM22.txt"   
mixture_file = 'LIHC_fpkm_mRNA_01A.txt'   #肿瘤患者表达谱
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
save(res_cibersort,file = "res_cibersort.Rdata")   #保存中间文件

load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22]   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细�?
#可视化（阿琛老师�?
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框写
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本�?
        yaxt = "n", # 先不绘制y�?
        ylab = "Relative percentage", # 修改y轴名�?
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板


####基因表达与ciber的相关�?####
setwd("cor")
install.packages("corrplot")
library(corrplot)
library(tidyverse)
expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- c("SPP1","PAGE1","G6PD","MAGEA4",'CDCA8',
          'TRIM54','KIF2C','KIF20A','ANLN',"SLC7A11")
exp <- expr[gene,]
exp <- exp %>% t() %>% as.data.frame()
ciber <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber <- ciber[,1:22]
identical(rownames(ciber),rownames(exp))
class(exp$SPP1)
class(ciber$`B cells naive`)
cor<-sapply(ciber,function(x,y) cor(x,y,method="spearman"),exp)
rownames(cor)<-colnames(exp)
cor_res <- cor.mtest(cor,#计算p�?
                     conf.level = 0.95)#置信区间
corrplot(cor,
         method = "color",#相关性矩阵展示的图形
         col=colorRampPalette(c("#01468b","white","#ee0000"))(100),
         addCoef.col = "black",#为相关系数添加颜�?
         tl.col="black",#设置文本标签的颜�?
         number.cex = 0.5,
         tl.cex = 0.7,
         cl.align = "l")
dev.off()


####根据基因高低组做生存分析####
setwd("survival")
surv <- read.table("exp_sur_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv$OS.time <- surv$OS.time*12
#median
median(surv$SPP1)
surv$group <- ifelse(surv$SPP1 > median(surv$SPP1),"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)

#3. 绘制生存曲线
#方法1
###3. 设置颜色，坐�?
plot(fit, conf.int = T,
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time(Months)",
     ylab = "Survival probablity(%)"
)
###添加标签
legend("topright",
       title = "Group",
       c("Low", "High"),
       lwd = 2, lty = 1,
       col = c("blue", "red"))
###添加P�?
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
text(25, 0.2, p.lab)
dev.off()


#方法2
install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险�?
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,120), # x轴长度，一般为0-10�?
           break.time.by = 20, # x轴步长为20个月
           legend.title = "",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标�?
           xlab = "Time (Months)", # 修改x轴标�?
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()

####SPP1差异分析####
setwd("SPP1_deg")
#install.packages("BiocManager") 
#if(!require(DESeq2))BiocManager::install('DESeq2')
library(DESeq2)
library(tidyverse)
counts_01A <- read.table("LIHC_counts_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- read.table("LIHC_fpkm_mRNA_01A.txt", sep = "\t",row.names = 1,check.names = F,header = T)
com <- intersect(colnames(counts_01A),colnames(exp))
exp <- exp[,com]
counts_01A <- counts_01A[,com]
identical(colnames(counts_01A),colnames(exp))
gene <- "SPP1"#每次运行只改这个基因�?
med=median(as.numeric(exp[gene,]))

conditions=data.frame(sample=colnames(exp),
                      group=factor(ifelse(exp[gene,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)

dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds)
save(res,file="res_deseq2_SPP1.Rda")
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)


####SPP1差异分析结果富集分析####
setwd("SPP1_fuji")
#install.packages("tidyverse")
#install.packages("BiocManager")
BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)

DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

#GO分析
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
save(ego,ego_res,file = "GO_SPP1_DEG.Rdata")

#3. 可视�?
##3.1 柱状�?
barplot(ego, showCategory = 20,color = "pvalue")
##3.2 气泡�?
dotplot(ego, showCategory = 20)
##3.3 分类展示
barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dotplot(ego,showCategory = 10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')


##### ConsensusClusterPlus聚类 #####
setwd("2_cluster")
#install.packages("BiocManager")
BiocManager::install('ConsensusClusterPlus')
library(tidyverse)
library(ConsensusClusterPlus)
exp <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
d=as.matrix(exp)
gene <- c("ALKBH5","YTHDF1") 
d <- d[gene,]
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:2],] #5是基因个�? 
d = sweep(d,1, apply(d,1,median,na.rm=T))

title=("JULEI") ##文件夹输出图片的位置
set.seed(1) #我发现设不设置种子都一�?
results = ConsensusClusterPlus(d,maxK=9,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1,plot="pdf")
results[[2]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]
results[[2]][["consensusClass"]][1:5]
icl = calcICL(results,title=title,plot="pdf") ##画另一组图�?

group<-results[[2]][["consensusClass"]]
group<-as.data.frame(group)
group$group <- factor(group$group,levels=c(1,2))
save(group,file = "group_AY.Rda")
load("group_GENE23456FINAL.Rda")

exp_gene <- exp[gene,]

# 绘制ConsensusClusterPlus后的热图
library(pheatmap)
group <- group %>% 
  rownames_to_column("sample")
annotation <- group %>% arrange(group) %>% column_to_rownames("sample")
a <- group %>% arrange(group) %>% mutate(sample=substring(.$sample,1,12))
b <- t(exp_gene) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  mutate(sample=substring(.$sample,1,12))
c <- inner_join(a,b,"sample") %>% .[,-2] %>% column_to_rownames("sample") %>% t(.)
pheatmap(c,annotation = annotation,
         cluster_cols = F,fontsize=5,fontsize_row=5,
         scale="row",show_colnames=F,
         fontsize_col=3)
pheatmap(c,annotation = annotation,
         annotation_colors = list(group = c("1" ="#01468b","2"= "#ee0000")),
         cluster_cols = F,fontsize=5,fontsize_row=5,
         scale="row",show_colnames=F,cluster_row = F,
         fontsize_col=3)
dev.off()

####cluster分组CIBERSORT画小提琴�?####
setwd("cibersort_AY")
library(tidyverse)
a <- read.table("CIBERSORT-Results.txt", sep = "\t",row.names = 1,check.names = F,header = T)
a <- a[,1:22]
identical(rownames(a),rownames(group))
b <- group
class(b$group)
a$group <- b$group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)
#install.packages("ggsci")
#install.packages("tidyr")
#install.packages("ggpubr")

b <- gather(a,key=CIBERSORT,value = Proportion,-c(group,sample))

ggboxplot(b, x = "CIBERSORT", y = "Proportion",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()

## ssGSEA
setwd("ssGSEA_AY")
BiocManager::install('GSVA')
library(tidyverse)
library(data.table)
library(GSVA)
#1.2 准备细胞marker
cellMarker <- data.table::fread("cellMarker.csv",data.table = F)
colnames(cellMarker)[2] <- "celltype"
#将cellMarker文件列名的第2个修改为celltype
type <- split(cellMarker,cellMarker$celltype)
#将cellMarker文件以celltype为分组拆分成list数据格式
#处理data.tables列表通常比使用group by参数按组对单个data.table进行操作要慢得多
cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})
#将list中每个celltype中的基因进行合并
save(cellMarker,file = "cellMarker_ssGSEA.Rdata")#保存中间文件
load("immune_infiltration//cellMarker_ssGSEA.Rdata")
##1.3 表达量矩阵的准备
###行是基因，列是样�?
expr <- data.table::fread("LIHC_fpkm_mRNA_01A.txt",data.table = F)   #读取表达文件
rownames(expr) <- expr[,1]   #将第一列作为行�?
expr <- expr[,-1]   #去除第一�?
expr <- as.matrix(expr)   #将expr转换为矩阵格�?

#2. 使用ssGSEA量化免疫浸润
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")

a <- gsva_data %>% t() %>% as.data.frame()
identical(rownames(a),rownames(group))
a$group <- group$group
a <- a %>% rownames_to_column("sample")
write.table(a,"ssGSEA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=ssGSEA,value = Expression,-c(group,sample))

ggboxplot(b, x = "ssGSEA", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()

####SPP1差异分析结果KEGG富集分析####
setwd("SPP1_KEGG")
#install.packages("tidyverse")
#install.packages("BiocManager")
BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)

DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

#KEGG分析
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk_res <- kk@result
save(kk,kk_res,file = "KEGG_SPP1_DEG.Rdata")

load("KEGG_SPP1_DEG.Rdata")

#柱状�?
barplot(kk, showCategory = 20,color = "pvalue")
#气泡�?
dotplot(kk, showCategory = 20)

dev.off()


####GSEA####
setwd("GSEA_SPP1")
#install.packages("tidyverse")
#install.packages("BiocManager")
BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)

DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))


msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "c5.all.v7.0.entrez.gmt"    #c2.all.v7.0.entrez.gmt �? c5.all.v7.0.entrez.gmt
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

geneList = DEG[,3]
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)

set.seed(1)
KEGG<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
#转换成数据框
KEGG_result_df <- as.data.frame(KEGG)
write.table(KEGG_result_df,file="GSEA_MSigDb_C5_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
save(KEGG,KEGG_result_df,file = "GSEA_deg_SPP1.rda")
#单个图绘�?
library(enrichplot)
gseaplot2(KEGG,1,color="red")
gseaplot2(KEGG,3,color="red",pvalue_table = T)

#汇总结�?
gseaplot2(KEGG, geneSetID = c(1,4,21,23,25,43), subplots = 1:3)
gseaplot2(KEGG, geneSetID = c(1,3), subplots = 1:3)
gseaplot2(KEGG, geneSetID = 1:3, subplots = 1)
gseaplot2(KEGG, geneSetID = 1:10, subplots = 1:3)

dev.off()


####LASSO COX回归模型的构�?####
setwd("lasso")

#install.packages("glmnet")
#install.packages("survival")
library("glmnet")
library("survival")

rt=read.table("data.exp.txt",header=T,sep="\t",row.names=1)           
rt$futime=rt$futime/365   


set.seed(3)   
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
plot(fit, xvar = "lambda", label = TRUE)

cvfit = cv.glmnet(x, y, family="cox", maxit = 1000)
plot(cvfit)
#其中两条虚线分别指示了两个特殊的λ�?
dev.off()
###4. 输出预测模型的相关系数与riskScore
###4.1 输出相关系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
geneCoef   #查看模型的相关系�?

###4.2 计算riskScore
FinalGeneExp = rt[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c("futime", "fustat", lassoGene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)

###5. 绘制散点分布�?
#install.packages("ggpubr")
library(ggpubr)  
p <- ggboxplot(dat, x = "fustat", y = "riskScore",
               color = "fustat", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   #得出预测结果

###6. 判断预测结果的准确�?
#install.packages("ROCR")
library(ROCR)   #使用ROCR包绘制预测模型的ROC曲线
library(glmnet)
library(caret)

pred <- prediction(dat$riskScore, dat$fustat)
perf <- performance(pred,"tpr","fpr")
AUC <- performance(pred,"auc")   #计算AUC
plot(perf,colorize=FALSE, col="red", print.auc =TRUE) #绘制ROC曲线
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
dev.off()

#画风险分布图
#y=生存时间
rt <- dat
color=as.vector(rt$fustat)
color[color==1]="indianred1"
color[color==0]="lightseagreen"
plot(rt$futime, pch=19,
     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("indianred1","lightseagreen"),cex=1.2)
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
abline(v=lowLength,lty=2)
dev.off()

#y=riskscore
rt <- rt[order(rt[,25]),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="Risk score",
     col=c(rep("lightseagreen",lowLength),rep("indianred1",highLength)) )
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("indianred1","lightseagreen"),cex=1.2)
dev.off()

####肿瘤突变负荷TMB####
#TCGA突变数据下载
#网址：https://portal.gdc.cancer.gov/
library(TCGAbiolinks)
#安装maftools
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("maftools")
library(maftools)
library(tidyverse)
setwd("TMB")
library(readxl)
library(readr)
#突变数据下载后，整理成数据框
mut2 <- read.maf("TCGA.LIHC.varscan.40fe9c1b-19d0-45cf-898a-f7b0cbad783e.DR-10.0.somatic.maf")

a <- mut2@data %>% 
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>% 
  as.data.frame() %>% 
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

gene <- as.character(unique(a$Hugo_Symbol))
sample <- as.character(unique(a$Tumor_Sample_Barcode))

mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))

for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))
gene_top <- gene_count$gene[1:20] # 修改数字，代表TOP多少
save(mat,mat_0_1,file = "TMB-LIHC.rda")

#以下绘图代码来自解螺旋阿琛老师
oncoplot(maf = mut2,
         top = 30,   #显示�?30个的突变基因信息
         fontSize = 0.6,   #设置字体大小
         showTumorSampleBarcodes = F)   #不显示病人信�?
dev.off()

####计算TMB####
maf = tmb(maf = mut2,
          captureSize = 50,
          logScale = TRUE)   
maf$sample <- substr(maf$Tumor_Sample_Barcode,1,16)
maf$sample <- gsub("-",".",maf$sample)
rownames(maf) <- maf$sample
#手动更改后导�?
write.csv(maf, file = "maf.csv")

com <- intersect(rownames(maf),rownames(group))
maf <- maf[com,]
group$x <- as.character(group$group)
group <- group %>% t() %>% as.data.frame()
group <- group[-1,]
group <- group[,com]
group <- group %>% t() %>% as.data.frame()
identical(rownames(group),rownames(maf))
tmb <- cbind(maf,group)
write.csv(tmb, file = "tmb.csv")


####nomogram####
setwd("nomogram")
library(tidyverse)
#整理数据
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
#读取表达数据
expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]
#双击差异分析文件
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)#根据自己需�?
#提取变化2倍的差异基因表达谱，进行后续lasso
expr <- expr[rownames(res_deseq2),]
expr <- expr %>% t() %>% as.data.frame()
#整合表达谱与生存信息
identical(rownames(expr),rownames(surv))
expr_surv <- cbind(surv,expr)
###lasso
#install.packages("glmnet")
#install.packages("survival")
library("glmnet")
library("survival")
colnames(expr_surv)[1] <- 'fustat'
colnames(expr_surv)[2] <- 'futime'
rt <- expr_surv         
rt$futime=rt$futime/365   

set.seed(3)   
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
plot(fit, xvar = "lambda", label = TRUE)

cvfit = cv.glmnet(x, y, family="cox", maxit = 1000)
plot(cvfit)
#其中两条虚线分别指示了两个特殊的λ�?
dev.off()
###4. 输出预测模型的相关系数与riskScore
###4.1 输出相关系数
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
geneCoef   #查看模型的相关系�?

###4.2 计算riskScore
FinalGeneExp = rt[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c("futime", "fustat", lassoGene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
dat = cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)

###5. 绘制散点分布�?
#install.packages("ggpubr")
library(ggpubr)  
p <- ggboxplot(dat, x = "fustat", y = "riskScore",
               color = "fustat", palette = "jco",
               add = "jitter")
p <- p + stat_compare_means()   #  Add p-value
p   #得出预测结果

###6. 判断预测结果的准确�?
#install.packages("ROCR")
library(ROCR)   #使用ROCR包绘制预测模型的ROC曲线
library(glmnet)
library(caret)

pred <- prediction(dat$riskScore, dat$fustat)
perf <- performance(pred,"tpr","fpr")
AUC <- performance(pred,"auc")   #计算AUC
plot(perf,colorize=FALSE, col="red", print.auc =TRUE) #绘制ROC曲线
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
dev.off()
#临床信息下载整理
#网址：https://xenabrowser.net/datapages/   下载后解压缩放入当前文件�?
#导入临床信息（不要行名）
clinical$submitter_id.samples <- gsub("-",".",clinical$submitter_id.samples)
rownames(clinical) <- clinical$submitter_id.samples
clinical <- clinical[,-1]
##结合dat与clinical
comgene1 <- intersect(rownames(clinical),rownames(dat))
dat <- dat[comgene1,]
clinical <- clinical[comgene1,]
dat <- dat[,-(3:13)]
clinical <- clinical[,c("gender.demographic","tumor_stage.diagnoses","pathologic_T","pathologic_N","pathologic_M")]
#clinical <- clinical1
identical(rownames(dat),rownames(clinical))
rt2 <- cbind(dat,clinical)
write.csv(rt2, file = "rt2.csv")
#导入修改后的rt3
rt2 <- rt3
#加载�?
library(rms)
library(foreign)
library(survival)
###3. 设置参数
rt2$gender <- factor(rt2$gender,labels=c("F", "M"))
rt2$stage <- factor(rt2$stage,labels=c("Stage1", "Stage2", "Stage3", "Stage4"))
rt2$T <- factor(rt2$T,labels=c("T1", "T2", "T3", "T4"))
rt2$M <- factor(rt2$M,labels=c("M0", "M1"))
rt2$N <- factor(rt2$N,labels=c("N0", "N1"))
rt2$risk <- factor(rt2$risk,labels=c("low", "high"))

ddist <- datadist(rt2)
options(datadist='ddist')   #使用函数datadist()将数据打�?

###4. 构建Cox回归模型
f <- cph(Surv(futime, fustat) ~gender + stage +T + M + N + risk, x=T, y=T, surv=T, data=rt2, time.inc=1)
surv <- Survival(f)

###5. 构建Nomogram
nom2 <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)), 
                 lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
                 maxscale=100, 
                 fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))
plot(nom2)
dev.off()


####合并GEO数据�?####
#网址：https://cloud.tencent.com/developer/article/1521695
setwd("batch")
#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("sva")
library(sva)
library(tidyverse)

load("GSE3325.Rda")
load("GSE46234.Rda")

merge_eset=inner_join(exprSet2,exprSet_GSE46234,by="symbol")
rownames(merge_eset) <- merge_eset$symbol
merge_eset <- merge_eset[,-1]
dim(merge_eset)
exp <- as.matrix(merge_eset)
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
dim(data)
class(data)
batchType <- c(rep(1,19),rep(2,8))
modType <- c(rep("normal",6),rep("tumor",13),rep("normal",4),rep("tumor",4))
mod  <-  model.matrix(~as.factor(modType))
outTab <- data.frame(ComBat(data, batchType,mod, par.prior=TRUE))

write.table(outTab,file="normalize.txt",sep="\t",quote=F,col.names=F)


####免疫检查点####
setwd("PD1")
library(tidyverse)
expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
a <- group %>% filter(group == "1")
b <- group%>% filter(group == "2")
expr1 <- expr[,rownames(a)]
expr2 <- expr[,rownames(b)]
gene <- c("PD1","PDL1","PDL2")
gene <- c("PDCD1","CD274","PDCD1LG2")
gene1 <- expr1[gene,]
gene2 <- expr2[gene,]
gene1 <- gene1 %>% t() %>% as.data.frame()
gene2 <- gene2 %>% t() %>% as.data.frame()
#https://www.genecards.org/
write.csv(gene1, file = "gene1.csv")
write.csv(gene2, file = "gene2.csv")


####fpkm转tpm####
setwd("TPM")
#定义转换函数
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
#读取
fpkm <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
fpkm <- 2^(fpkm)-1
tpm <- as.data.frame(round(apply(fpkm,2,fpkmToTpm),2))

####X cell####
#安装devtools 右下角install
##安装Rtools�?
#b站教�? 视频号：BV1bk4y1k7Nd 请认真看�? 给这位up一键三�?
#下载网址windows版本 https://cran.r-project.org/bin/windows/Rtools/
#验证Rtools安装是否成功
system("g++ -v")
system("where make")
#安装不成功请仔细看b站教�? 视频号：BV1bk4y1k7Nd
#安装xCell�?
chooseBioCmirror()
devtools::install_github('dviraran/xCell')
library(xCell)
library(ggpubr)
library(tidyverse)
setwd("Xcell")
exp <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#查看细胞类型
celltypeuse<-xCell.data$spill$K
rs<-xCellAnalysis(exp,parallel.sz=10) #计算
#准备分组信息
gene <- "PDCD1"#每次运行只改这个基因�?
med=median(as.numeric(exp[gene,]))
expgene <- exp[gene,]
expgene <- expgene %>% t() %>% as.data.frame()
expgene$group <- ifelse(expgene$PDCD1>med,"High","Low") #�?
rs <- as.data.frame(rs)
rs <- rs %>% t() %>% as.data.frame()
comname <- intersect(rownames(rs),rownames(expgene)) 
rs <- rs[comname,]
expgene <- expgene[comname,]
identical(rownames(rs),rownames(expgene))
rs$group <- as.factor(expgene$group)
class(rs$group)
rs <- rs %>% rownames_to_column("sample")
a <- rs
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=xCell,value = Expression,-c(group,sample))

ggboxplot(b, x = "xCell", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()

