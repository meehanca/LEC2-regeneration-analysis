library(DESeq2)
library(ggplot2)
library(gplots)
library(reshape2)
library(pheatmap)
library(VennDiagram)
library(ggrepel)


setwd("~/Documents/WGCNA/LEC2/transcriptomics/LEC2 counts")
sampleDataFilename <- 'sampleTable_renamed.txt'

sampleTable = read.table(sampleDataFilename,header=TRUE)

head(sampleTable)
htseqDir<-getwd()

##  Read in the results from the LibiNorm analysis (the counts files)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = htseqDir,design = ~  genotype)

##  And perform the analysis (details in the manual)
dds<-DESeq(ddsHTSeq)

####################################################################
# Do PCA
####################################################################
#principal component analysis

vst = vst(dds)

v <- plotPCA(vst, intgroup=c("genotype"))
v<- v+ geom_label_repel(aes(label = name))
v
pcaData <- DESeq2::plotPCA(vst, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA_parents.pdf", height = 6, width = 6)
ggplot(pcaData, aes(PC1, PC2, color=genotype, shape=genotype)) +
  geom_point(size=3) +
  #scale_colour_manual(name="",values = c("a12"="goldenrod2", "gd33"="darkslateblue", "f1"="saddlebrown"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_bw()
#dev.off()

####################################################################
#Plotting Reps
####################################################################
plot_reps =  function(dds,x=1,y=2,cond_choice=1, cond='genotype'){
  ##  Estimate the size factors for normalisation
  dds<-estimateSizeFactors(dds)
  
  ## Extract the normalised counts for the condition you want
  rep_values<- counts(dds, normalized=TRUE)[,dds[[cond]]==cond_choice]
  
  # Take logs of these values
  vals <- log2(rep_values[,c(x,y)] + 0.5)
  # And plot
  plot(vals,pch=16, cex=0.4,xlab=paste('rep',x),ylab=paste('rep',y))
  grid(col = "darkgray", lty = "solid",lwd = par("lwd"), equilogs = TRUE)
  title(paste("Comparison of",cond_choice,"replicates"))
}

####################################################################
#DEGs
####################################################################
filter_degs <- function(res){
  summary(res)
  res2 = res[!(is.na(res$padj)),]
  res2 = res2[res2$padj < 0.01,]
  return(res2)
}

resultsNames(dds)
LO_root_DEGs = results(dds, contrast= c("genotype", "LO_root", "LEC2_root"), alpha = 0.01, pAdjustMethod = "BH", lfcThreshold = 0)
LO_root_DEG = filter_degs(LO_root_DEGs)

LO_leaf_DEGs = results(dds, contrast= c("genotype", "LO_leaf", "LEC2_leaf"), alpha = 0.01, pAdjustMethod = "BH", lfcThreshold = 0)
LO_leaf_DEG = filter_degs(LO_leaf_DEGs)

RO_root_DEGs = results(dds, contrast= c("genotype", "RO_root", "LEC2_root"), alpha = 0.01, pAdjustMethod = "BH", lfcThreshold = 0)
RO_root_DEG= filter_degs(RO_root_DEGs)

RO_leaf_DEGs = results(dds, contrast= c("genotype", "RO_leaf", "LEC2_leaf"), alpha = 0.01, pAdjustMethod = "BH", lfcThreshold = 0)
RO_leaf_DEG = filter_degs(RO_leaf_DEGs)

###################################################################
#ATG WS-2 -> Col-0
###################################################################

WS_2_Col_0 <- read.delim('WS-2_Col-0.txt')

rownames(LO_root_DEG) <- as.vector(merge(rownames(LO_root_DEG),
                                         WS_2_Col_0, by.x=1, by.y=1)[,2])

rownames(LO_leaf_DEG) <- as.vector(merge(rownames(LO_leaf_DEG),
                                         WS_2_Col_0, by.x=1, by.y=1)[,2])

rownames(RO_leaf_DEG) <- as.vector(merge(rownames(RO_leaf_DEG),
                                         WS_2_Col_0, by.x=1, by.y=1)[,2])

rownames(RO_root_DEG) <- as.vector(merge(rownames(RO_root_DEG),
                                         WS_2_Col_0, by.x=1, by.y=1)[,2])

###################################################################
#Preparing lists of diffrentially expressed genes for online tools
####################################################################

LO_root_up_DEG <- LO_root_DEG[LO_root_DEG$log2FoldChange>0,]
LO_root_down_DEG <- LO_root_DEG[LO_root_DEG$log2FoldChange<0,]

LO_leaf_up_DEG <- LO_leaf_DEG[LO_leaf_DEG$log2FoldChange>0,]
LO_leaf_down_DEG <- LO_leaf_DEG[LO_leaf_DEG$log2FoldChange<0,]

RO_root_up_DEG <- RO_root_DEG[RO_root_DEG$log2FoldChange>0,]
RO_root_down_DEG <- RO_root_DEG[RO_root_DEG$log2FoldChange<0,]

RO_leaf_up_DEG <- RO_leaf_DEG[RO_leaf_DEG$log2FoldChange>0,]
RO_leaf_down_DEG <- RO_leaf_DEG[RO_leaf_DEG$log2FoldChange<0,]

###################################################################
#Preparing lists of diffrentially expressed genes for online tools
####################################################################

write.table(rownames(LO_root_up_DEG), file= 'LO_root_up_DEG.txt', quote =F, row.names = F, col.names = F)
write.table(rownames(LO_root_down_DEG), file= 'LO_root_down_DEG.txt', quote =F, row.names = F, col.names = F)

write.table(rownames(LO_leaf_up_DEG), file= 'LO_leaf_up_DEG.txt', quote =F, row.names = F, col.names = F)
write.table(rownames(LO_leaf_down_DEG), file= 'LO_leaf_down_DEG.txt', quote =F, row.names = F, col.names = F)

write.table(rownames(RO_root_up_DEG), file= 'RO_root_up_DEG.txt', quote =F, row.names = F, col.names = F)
write.table(rownames(RO_root_down_DEG), file= 'RO_root_down_DEG.txt', quote =F, row.names = F, col.names = F)

write.table(rownames(RO_leaf_up_DEG), file= 'RO_leaf_up_DEG.txt', quote =F, row.names = F, col.names = F)
write.table(rownames(RO_leaf_down_DEG), file= 'RO_leaf_down_DEG.txt', quote =F, row.names = F, col.names = F)
