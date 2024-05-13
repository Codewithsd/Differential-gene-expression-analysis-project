####..Diffrential Gene Expression Analysis for papillary Thyroid carcinoma ..####

#load required libraries 
if (!requireNamespace("zoo", quietly = TRUE)) 
  install.packages("zoo")
if(!require("readxl",quietly = TRUE))
  install.packages("readxl")
if(!require("ggplot",quietly=TRUE))
  install.packages("ggplot",force=TRUE)
BiocManager::install("apeglm",force = TRUE)

library(apeglm)
library(DESeq2)
library(tidyverse)
library(zoo)
library(pheatmap)
library(EnhancedVolcano)



#load data
#load raw count data 

counts_data<-read.delim('C:/Users/santo/OneDrive/Desktop/project/RNA-Sequencing of human papillary thyroid carcinomas/E-GEOD-64912-raw-counts 2.tsv',row.names = 1)
dim(counts_data)
head(counts_data)


#load sample info table

sample_info<- read.table('C:/Users/santo/OneDrive/Desktop/project/RNA-Sequencing of human papillary thyroid carcinomas/sample info.txt',sep="\t",header =TRUE,dec = ".")
dim(sample_info)
head(sample_info)


#check the column of counts table are same as rownames in in sample_info 

all(rownames(sample_info) == colnames(counts_data))


#data pre processing 


counts_data<- subset(counts_data,select=-Gene.Name )


#preparing sample info table 

rownames(sample_info)<- sample_info$Run
sample_info$Run<-NULL

dim(sample_info)

dim(counts_data)

all(rownames(sample_info) == colnames(counts_data))

counts_filterd<-counts_data[rowSums(counts_data)>=10,]

#set factor levels
factors<-factor(sample_info$type)
groups <- unique(sample_info$type) 
groups

groups<-rev(groups)

sample_info$type<-factors

sample_info$type



#########   DEseq2 analysis   ########## 


#create Deseq Object
dds<-DESeqDataSetFromMatrix(countData = counts_filterd,
                            colData = sample_info,
                            design = ~type)

dds$type<-relevel(dds$type,ref = 'normal')
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# Extract shrinkage dispersion estimates
disp <- dispersionFunction(dds)

# Optionally, you can plot dispersion estimates
plotDispEsts(dds)

#filter out low read counts <10
#keep<- rowSums(counts(dds)>=10)<=min(table(sample_info$type))
#dds<-dds[keep,]
#head(keep)

#performing deseq function 
dds<-DESeq(dds,test="Wald",sfType='poscount')
head(dds)



#save normalised counts from dds object matrix 
normalized_counts<-counts(dds,normalized=TRUE)
head(normalized_counts)
write.csv(normalized_counts,"normalised_counts.csv")

res<-results(object = dds, contrast = c('type','papillary_thyroid_carcinoma','normal'),
                      pAdjustMethod = 'BH',alpha = 0.05)
summary(res)

res <- results(dds, contrast=c('type','papillary_thyroid_carcinoma','normal'))
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
summary(res)

#res<-results(dds)
#making result as a data frame

res 

res<-as.data.frame(res)
class(res)
head(res) 

dim(res)
names(res)

#res <- results(dds, name="normal vs papillary_thyroid_carcinoma")
#res <- results(dds, contrast=c('type','papillary_thyroid_carcinoma','normal))
#resLFC <- lfcShrink(dds, coef="normal vs papillary_thyroid_carcinoma", type="apeglm")


#introdusing gene id as a new column
res$Gene.ID<-row.names(res)
names(res) 
head(res)


res<-subset(res,
            select=c("Gene.ID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj" )) 



reorderd<- res[order(res$padj),]
reorderd
#extract de genes with padje<0.05 and log2foldchange <=-1 or >=1

deg<-subset(res,padj<0.05) # & abs(log2FoldChange)>= 1)

dim(deg)

head(deg)

deg<- deg[order(deg$padj),]
deg<-as.data.frame(deg)
novalgene<- select(deg,c(Gene.ID))
write.csv(novalgene,"novel_genes.csv")

#Variance stabilizing transformation

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

###########   VISUALISATIONS  ############

#Dispersion plot for dispersion estimate
plotDispEsts(dds)

#histogram plot of p-values
hist(res$padj,breaks = seq(0,1,length=20, col="grey",border="white",
                           xlab="",ylab="",ylim=c(0,6000),main="histogram for p values "))


#Maplot 
limma::plotMA(res,cex=0.7,ylim=c(-10,10))
abline(h=c(-1,1,col="red",lwd=1))
plotMA(res,alpha = 0.05)
abline(h=c(-2,2),col = "blue")




#volcano plot for diffrentialy expressed genes 
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')


#heatmap
res <- as.data.frame(res)
top_genes <- res[which(res$padj < 0.05), ]
heatmap_data <- counts[rownames(top_genes), ]
heatmap(heatmap_data)

#pca plot

plotPCA(rld,intgroup=c("type"))



##### GSEW #####


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler",force = TRUE)
library(clusterProfiler)

BiocManager::install("AnnotationDbi",force = TRUE)

BiocManager::install("org.Hs.eg.db")

library(AnnotationDbi)
library(org.Hs.eg.db)

#preparing genes for enrichment analysis 
genes<- rownames(novalgene)
write.table(genes,"genes.csv",sep=",")

#using enrichGO function 

GO_result<- enrichGO(gene=genes,OrgDb = "org.Hs.eg.db", keyType="ENSEMBL", ont="BP")

as.data.frame(GO_result)

fit<- plot(barplot(GO_result,showCatagory=20))

fit

library(enrichplot)
goplot(ego)

