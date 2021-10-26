#http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#preparing-count-matrices consulted 07.03.2016

#source("https://bioconductor.org/biocLite.R")
#biocLite(c("EDASeq","DESeq2","pheatmap","ggplot2","gplots","geneplotter","RColorBrewer","genefilter"))

library("EDASeq")
library("DESeq2")
library("pheatmap")
library("ggplot2")
library("gplots")
library("geneplotter")
library("RColorBrewer")
library("genefilter")

library("BiocParallel")
register(MulticoreParam(4))

library("dplyr")


#Import information of experimental setup

CORAL_INFO<-read.csv("PATH TO SAMPLE INFO FILE", head=T)

# Import counts matrix (metatranscriptome)

metatranscriptome_count_table <- read.table("PATH TO FILE", header = TRUE, sep = '\t')


#Import host sequences IDs (determined with psytrans)

host_sequences_ids <- read.table('PATH TO FILE', header=FALSE, sep = '\t')


#Obtain the host counts from the metatranscriptome table

CORAL_DE <- subset(metatranscriptome_count_table, rownames(metatranscriptome_count_table) %in% host_sequences_ids$V1)


#Create a 'not in' operator to extract all sequences NOT present in host IDs list

`%notin%` <- Negate(`%in%`)


#Obtain the symbiont counts from the metatranscriptome table.

Symbiont_DE <- subset(metatranscriptome_count_table, rownames(metatranscriptome_count_table) %notin% host_sequences_ids$V1)


					

# Check columns names 


colnames(HOST_DE)
colnames(Symbiont_DE)


# Remove genes that have no counts over all samples
HOST_DE <- HOST_DE[(rowSums(HOST_DE) > 0),]
Symbiont_DE <- Symbiont_DE[(rowSums(Symbiont_DE) > 0),]

HOST_DE <- floor(HOST_DE)
HOST_DESeq<-DESeqDataSetFromMatrix(countData=HOST_DE, colData=HOST_DE_INFO, design=~experiment + condition)

Symbiont_DE <- floor(Symbiont_DE)
Symbiont_DESeq<-DESeqDataSetFromMatrix(countData=Symbiont_DE, colData=CORAL_INFO, design=~experiment + condition)


#relevel the factors so that control is the reference
HOST_DESeq$condition<-relevel(HOST_DESeq$condition,ref="control")
Symbiont_DESeq$condition<-relevel(Symbiont_DESeq$condition,ref="control")

#estimate size factors
#Thus, if all size factors are roughly equal to one, the libraries have been sequenced equally deeply.
HOST_DESeq<-estimateSizeFactors(HOST_DESeq)
Symbiont_DESeq<-estimateSizeFactors(Symbiont_DESeq)

#get rows with all non-zero counts
non_zero_rows<-apply(counts(HOST_DESeq), 1, function(x){all(x>0)})
non_zero_rows<-apply(counts(Symbiont_DESeq), 1, function(x){all(x>0)})

#get all rows
all_rows<-apply(counts(HOST_DESeq), 1, function(x){all(x>=0)})
all_rows<-apply(counts(Symbiont_DESeq), 1, function(x){all(x>=0)})

#number of non-zero rows:
sum(non_zero_rows)

#number of rows:
sum(all_rows)

#cummulative distribution of normalized counts for non-zero rows
multiecdf(counts(HOST_DESeq, normalized=T)[non_zero_rows, ], xlab="mean counts", xlim=c(0,1000))
multiecdf(counts(Symbiont_DESeq, normalized=T)[non_zero_rows, ], xlab="mean counts", xlim=c(0,1000))

#density of normalized counts
multidensity(counts(HOST_DESeq, normalized=T)[non_zero_rows, ], xlab="mean counts", xlim=c(0,1000))
multidensity(counts(Symbiont_DESeq, normalized=T)[non_zero_rows, ], xlab="mean counts", xlim=c(0,1000))

#sample comparisons
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(1,2), main=paste( colnames(HOST_DESeq)[1], "vs.", colnames(HOST_DESeq)[2]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(1,3), main=paste( colnames(HOST_DESeq)[1], "vs.", colnames(HOST_DESeq)[3]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(1,4), main=paste( colnames(HOST_DESeq)[1], "vs.", colnames(HOST_DESeq)[4]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(1,5), main=paste( colnames(HOST_DESeq)[1], "vs.", colnames(HOST_DESeq)[5]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(1,6), main=paste( colnames(HOST_DESeq)[1], "vs.", colnames(HOST_DESeq)[6]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(1,7), main=paste( colnames(HOST_DESeq)[1], "vs.", colnames(HOST_DESeq)[7]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(1,8), main=paste( colnames(HOST_DESeq)[1], "vs.", colnames(HOST_DESeq)[8]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(2,3), main=paste( colnames(HOST_DESeq)[2], "vs.", colnames(HOST_DESeq)[3]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(2,4), main=paste( colnames(HOST_DESeq)[2], "vs.", colnames(HOST_DESeq)[4]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(2,5), main=paste( colnames(HOST_DESeq)[2], "vs.", colnames(HOST_DESeq)[5]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(2,6), main=paste( colnames(HOST_DESeq)[2], "vs.", colnames(HOST_DESeq)[6]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(2,7), main=paste( colnames(HOST_DESeq)[2], "vs.", colnames(HOST_DESeq)[7]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(2,8), main=paste( colnames(HOST_DESeq)[2], "vs.", colnames(HOST_DESeq)[8]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(2,9), main=paste( colnames(HOST_DESeq)[2], "vs.", colnames(HOST_DESeq)[9]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(3,4), main=paste( colnames(HOST_DESeq)[3], "vs.", colnames(HOST_DESeq)[4]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(3,5), main=paste( colnames(HOST_DESeq)[3], "vs.", colnames(HOST_DESeq)[5]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(3,6), main=paste( colnames(HOST_DESeq)[3], "vs.", colnames(HOST_DESeq)[6]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(3,7), main=paste( colnames(HOST_DESeq)[3], "vs.", colnames(HOST_DESeq)[7]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(3,8), main=paste( colnames(HOST_DESeq)[3], "vs.", colnames(HOST_DESeq)[8]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(3,9), main=paste( colnames(HOST_DESeq)[3], "vs.", colnames(HOST_DESeq)[9]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(4,5), main=paste( colnames(HOST_DESeq)[4], "vs.", colnames(HOST_DESeq)[5]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(4,6), main=paste( colnames(HOST_DESeq)[4], "vs.", colnames(HOST_DESeq)[6]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(4,7), main=paste( colnames(HOST_DESeq)[4], "vs.", colnames(HOST_DESeq)[7]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(4,8), main=paste( colnames(HOST_DESeq)[4], "vs.", colnames(HOST_DESeq)[8]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(4,9), main=paste( colnames(HOST_DESeq)[4], "vs.", colnames(HOST_DESeq)[9]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(5,6), main=paste( colnames(HOST_DESeq)[5], "vs.", colnames(HOST_DESeq)[6]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(5,7), main=paste( colnames(HOST_DESeq)[5], "vs.", colnames(HOST_DESeq)[7]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(5,8), main=paste( colnames(HOST_DESeq)[5], "vs.", colnames(HOST_DESeq)[8]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(5,9), main=paste( colnames(HOST_DESeq)[5], "vs.", colnames(HOST_DESeq)[9]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(6,7), main=paste( colnames(HOST_DESeq)[6], "vs.", colnames(HOST_DESeq)[7]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(6,8), main=paste( colnames(HOST_DESeq)[6], "vs.", colnames(HOST_DESeq)[8]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(6,9), main=paste( colnames(HOST_DESeq)[6], "vs.", colnames(HOST_DESeq)[9]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(7,8), main=paste( colnames(HOST_DESeq)[7], "vs.", colnames(HOST_DESeq)[8]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(7,9), main=paste( colnames(HOST_DESeq)[7], "vs.", colnames(HOST_DESeq)[9]))
#MDPlot(counts(HOST_DESeq, normalized=T)[non_zero_rows, ],c(8,9), main=paste( colnames(HOST_DESeq)[8], "vs.", colnames(HOST_DESeq)[9]))


#transform the data to rlog
HOST_DESeq_rlog<-rlogTransformation(HOST_DESeq, blind = T)
Symbiont_DESeq_rlog<-rlogTransformation(Symbiont_DESeq, blind = T)


#produce a heat plot using the transformed distances
HOST_DESeq_distances<-dist(t(assay(HOST_DESeq_rlog)))
Symbiont_DESeq_distances<-dist(t(assay(Symbiont_DESeq_rlog)))

heatmap.2(as.matrix(HOST_DESeq_distances), trace="none", col=rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))
heatmap.2(as.matrix(Symbiont_DESeq_distances), trace="none", col=rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))

#PCA of samples
plotPCA(HOST_DESeq_rlog, intgroup=c("condition"))
plotPCA(Symbiont_DESeq_rlog, intgroup=c("condition"))

plotPCA(HOST_DESeq_rlog, intgroup=c("experiment"))
plotPCA(Symbiont_DESeq_rlog, intgroup=c("experiment"))

#remove outlier samples: e.g. Cauris_141021_6T
outliers=c("")
HOST_DESeq_NoOutliers<-HOST_DESeq[,!(colnames(HOST_DESeq) %in% outliers)]
Symbiont_DESeq_NoOutliers<-Symbiont_DESeq[,!(colnames(Symbiont_DESeq) %in% outliers)]

#repeat above analyses
HOST_DESeq_NoOutliers_rlog<-rlogTransformation(HOST_DESeq_NoOutliers, blind = T)
Symbiont_DESeq_NoOutliers_rlog<-rlogTransformation(Symbiont_DESeq_NoOutliers, blind = T)


#produce a heat plot using the transformed distances
HOST_distances_no_outliers<-dist(t(assay(HOST_DESeq_NoOutliers_rlog)))
Symbiont_distances_no_outliers<-dist(t(assay(Symbiont_DESeq_NoOutliers_rlog)))

heatmap.2(as.matrix(HOST_distances_no_outliers), trace="none", col=rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))
heatmap.2(as.matrix(Symbiont_distances_no_outliers), trace="none", col=rev(colorRampPalette(brewer.pal(9, "GnBu"))(100)))

#PCA of samples
plotPCA(HOST_DESeq_NoOutliers_rlog, intgroup=c("condition"))
plotPCA(HOST_DESeq_NoOutliers_rlog, intgroup=c("sample_name"))
plotPCA(HOST_DESeq_NoOutliers_rlog, intgroup=c("experiment"))

#heat map of most variable genes
HOST_most_variable_genes<-head(order(rowVars(assay(HOST_DESeq_NoOutliers_rlog)), decreasing = T), 750)
Symbiont_most_variable_genes<-head(order(rowVars(assay(Symbiont_DESeq_NoOutliers_rlog)), decreasing = T), 750)

heatmap.2(assay(HOST_DESeq_NoOutliers_rlog)[HOST_most_variable_genes, ], scale="row", trace = "none", dendrogram = "column", col= colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

heatmap.2(assay(Symbiont_DESeq_NoOutliers_rlog)[Symbiont_most_variable_genes, ], scale="row", trace = "none", dendrogram = "column", col= colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))


###
#
#Differential expression analysis
#
#####
HOST_DESeq_NoOutliers<-estimateSizeFactors(HOST_DESeq_NoOutliers)
Symbiont_DESeq_NoOutliers<-estimateSizeFactors(Symbiont_DESeq_NoOutliers)

HOST_DESeq_NoOutliers<-estimateDispersions(HOST_DESeq_NoOutliers)
Symbiont_DESeq_NoOutliers<-estimateDispersions(Symbiont_DESeq_NoOutliers)

plotDispEsts(HOST_DESeq_NoOutliers)
plotDispEsts(Symbiont_DESeq_NoOutliers)

#wald test
HOST_DESeq_NoOutliers<-nbinomWaldTest(HOST_DESeq_NoOutliers)
HOST_DESeq_NoOutliers_Results<-results(HOST_DESeq_NoOutliers, pAdjustMethod = "BH", contrast=c('condition','control','TREATMENT'))

Symbiont_DESeq_NoOutliers<-nbinomWaldTest(Symbiont_DESeq_NoOutliers)
Symbiont_DESeq_NoOutliers_Results<-results(Symbiont_DESeq_NoOutliers, pAdjustMethod = "BH", contrast=c('condition','control','TREATMENT'))


#number of DE genes at 0.01 significance level
table(HOST_DESeq_NoOutliers_Results$padj < 0.01)
table(Symbiont_DESeq_NoOutliers_Results$padj < 0.01)

#number of DE genes at 0.01 significance level and log fold change > 2
table(HOST_DESeq_NoOutliers_Results$padj < 0.01 & abs(HOST_DESeq_NoOutliers_Results$log2FoldChange) > 2)
table(Symbiont_DESeq_NoOutliers_Results$padj < 0.01 & abs(Symbiont_DESeq_NoOutliers_Results$log2FoldChange) > 2)

#number of overexpressed genes at 0.01 significance level and log fold change > 2
table(HOST_DESeq_NoOutliers_Results$padj < 0.01 & HOST_DESeq_NoOutliers_Results$log2FoldChange > 2)
table(Symbiont_DESeq_NoOutliers_Results$padj < 0.01 & Symbiont_DESeq_NoOutliers_Results$log2FoldChange > 2)

#number of underexpressed genes at 0.01 significance level and log fold change < -2
table(HOST_DESeq_NoOutliers_Results$padj < 0.01 & HOST_DESeq_NoOutliers_Results$log2FoldChange < -2)
table(Symbiont_DESeq_NoOutliers_Results$padj < 0.01 & Symbiont_DESeq_NoOutliers_Results$log2FoldChange < -2)


#write results matrix to file
write.csv(HOST_DESeq_NoOutliers_Results, "PATH TO SAVE FILE",quotes=FALSE)
write.csv(Symbiont_DESeq_NoOutliers_Results, "PATH TO SAVE FILE",quotes=FALSE)


plotMA(HOST_DESeq_NoOutliers_Results, alpha=0.01)
plotMA(Symbiont_DESeq_NoOutliers_Results, alpha=0.01)

#get significant genes
degs_names_host <-rownames(subset(HOST_DESeq_NoOutliers_Results, padj<0.01 & abs(log2FoldChange) > 2))
write.csv(degs_names_host, "PATH TO SAVE FILE")

degs_names_symbiont <-rownames(subset(Symbiont_DESeq_NoOutliers_Results, padj<0.01 & abs(log2FoldChange) > 2))
write.csv(degs_names_host, "PATH TO SAVE FILE")

#get significantly overexpressed in treatment:
over_deg_names_host<-rownames(subset(HOST_DESeq_NoOutliers_Results, padj<0.01 & log2FoldChange > 2))
over_deg_names_symbiont<-rownames(subset(Symbiont_DESeq_NoOutliers_Results, padj<0.01 & log2FoldChange > 2))

write.csv(over_deg_names_host, "PATH TO SAVE FILE",quote=FALSE, row.names=FALSE)
write.csv(over_deg_names_symbiont, "PATH TO SAVE FILE",quote=FALSE, row.names=FALSE)

heatmap.2(log2(counts(HOST_DESeq_NoOutliers, normalized=T)[rownames(HOST_DESeq_NoOutliers) %in% over_deg_names_host, ]+1), scale="row", trace = "none",dendrogram = "column", col= colorRampPalette(rev(brewer.pal(9, "RdBu"))))

heatmap.2(log2(counts(Symbiont_DESeq_NoOutliers, normalized=T)[rownames(Symbiont_DESeq_NoOutliers) %in% over_deg_names, ]+1), scale="row", trace = "none", dendrogram = "column", col= colorRampPalette(rev(brewer.pal(9, "RdBu"))))

#get significantly underexpressed in treatment:
under_deg_names_host <-rownames(subset(HOST_DESeq_NoOutliers_Results, padj<0.01 & log2FoldChange < -2))
under_deg_names_symbiont<-rownames(subset(Symbiont_DESeq_NoOutliers_Results, padj<0.01 & log2FoldChange < -2))

write.csv(under_deg_names_host, "PATH TO SAVE FILE",row.names=FALSE, quote=FALSE)
write.csv(under_deg_names_symbiont, "PATH TO SAVE FILE",row.names=FALSE, quote=FALSE)

heatmap.2(log2(counts(HOST_DESeq_NoOutliers, normalized=T)[rownames(HOST_DESeq_NoOutliers) %in% under_deg_names, ]+1), scale="row", trace = "none", dendrogram = "column", col= colorRampPalette(rev(brewer.pal(9, "RdBu"))))

heatmap.2(log2(counts(Symbiont_DESeq_NoOutliers, normalized=T)[rownames(Symbiont_DESeq_NoOutliers) %in% under_deg_names, ]+1), scale="row", trace = "none", dendrogram = "column", col= colorRampPalette(rev(brewer.pal(9, "RdBu"))))
