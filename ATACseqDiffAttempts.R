# we could not get our code to work to run Differential analysis on the ATAC Seq data

setwd("/class/projects/bioinf545-2021/group-1/ATACseq")

library(Rsamtools)
library(ggplot2)
library(magrittr)
library(Rsubread)

indexBam("SRR13169286.pruned.bam")
indexBam("SRR13169287.pruned.bam")
indexBam("SRR13169288.pruned.bam")
indexBam("SRR13169289.pruned.bam")

sortedBAM286 <- "SRR13169286.pruned.bam"
sortedBAM287 <- "SRR13169287.pruned.bam"
sortedBAM288 <- "SRR13169288.pruned.bam"
sortedBAM289 <- "SRR13169289.pruned.bam"


 
#plot mapped reads
barplot(idxstatsBam(sortedBAM286)$mapped, horiz = T, ylab = "Chromosome", xlab = 'Number of Mapped Reads',
        names.arg = idxstatsBam(sortedBAM286)$seqnames, 
        main = "Stimulated ATAC seq data SRR13169286",  xlim = c(0, 2e+06))

barplot(idxstatsBam(sortedBAM287)$mapped, horiz = T, ylab = "Chromosome", 
        xlab = 'Number of Mapped Reads',
        names.arg = idxstatsBam(sortedBAM287)$seqnames, 
        main = "Stimulated ATAC seq data SRR13169287",  xlim = c(0, 2e+06))

barplot(idxstatsBam(sortedBAM288)$mapped, horiz = T, ylab = "Chromosome", 
        xlab = 'Number of Mapped Reads',
        names.arg = idxstatsBam(sortedBAM288)$seqnames, 
        main = "Stimulated ATAC seq data SRR13169288",  xlim = c(0, 2e+06))

barplot(idxstatsBam(sortedBAM289)$mapped, horiz = T, ylab = "Chromosome", 
        xlab = 'Number of Mapped Reads',
        names.arg = idxstatsBam(sortedBAM289)$seqnames, 
        main = "Stimulated ATAC seq data SRR13169289",  xlim = c(0, 2e+06))

#reading mapped reads

library(GenomicAlignments)

atacReads286 <- readGAlignmentPairs(sortedBAM286, param= ScanBamParam(mapqFilter = 1, 
                                                                      what = c("qname","mapq", "isize"),
                                                                      which = GRanges("1", IRanges(1,6302552))))
                                    
atacReads286_read1 <- GenomicAlignments::first(atacReads286)
insertSizes286 <- abs(elementMetadata(atacReads286_read1)$isize)



atacReads287 <- readGAlignmentPairs(sortedBAM287, param= ScanBamParam(mapqFilter = 1, 
                                                                      what = c("qname","mapq", "isize"),
                                                                      which = GRanges("1", IRanges(1,6302552))))

atacReads287_read1 <- GenomicAlignments::first(atacReads287)
insertSizes287 <- abs(elementMetadata(atacReads287_read1)$isize)


atacReads288 <- readGAlignmentPairs(sortedBAM288, param= ScanBamParam(mapqFilter = 1, 
                                                                      what = c("qname","mapq", "isize"),
                                                                      which = GRanges("1", IRanges(1,6302552))))

atacReads288_read1 <- GenomicAlignments::first(atacReads288)
insertSizes288 <- abs(elementMetadata(atacReads288_read1)$isize)

atacReads289 <- readGAlignmentPairs(sortedBAM289, param= ScanBamParam(mapqFilter = 1, 
                                                                      what = c("qname","mapq", "isize"),
                                                                      which = GRanges("1", IRanges(1,6302552))))

atacReads289_read1 <- GenomicAlignments::first(atacReads289)
insertSizes289 <- abs(elementMetadata(atacReads289_read1)$isize)

library(magrittr)
library(dplyr)
library(ggplot2)

table(insertSizes286) %>% data.frame %>% rename(InsertSize = insertSizes286, Count = Freq) %>% 
  mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count))) %>% 
  plot(aes(x= InsertSize, y = Count), type ="l")

table(insertSizes287) %>% data.frame %>% rename(InsertSize = insertSizes287, Count = Freq) %>% 
  mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count))) %>% 
  plot(aes(x= InsertSize, y = Count), type ="l")

table(insertSizes288) %>% data.frame %>% rename(InsertSize = insertSizes288, Count = Freq) %>% 
  mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count))) %>% 
  plot(aes(x= InsertSize, y = Count), type ="l")

table(insertSizes289) %>% data.frame %>% rename(InsertSize = insertSizes289, Count = Freq) %>% 
  mutate(InsertSize = as.numeric(as.vector(InsertSize)), Count = as.numeric(as.vector(Count))) %>% 
  plot(aes(x= InsertSize, y = Count), type ="l")

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(soGGi)
TSSs <- resize(genes(TxDb.Mmusculus.UCSC.mm10.knownGene),fix = "start",1)

nucFree286 <- regionPlot(bamFile = sortedBAM286, testRanges = TSSs, style = "point", 
                      format = "bam", paired = TRUE, minFragmentLength = 0, maxFragmentLength = 100, 
                      forceFragment = 50)

atacReads286_Open <- atacReads286[insertSizes286 < 100,]
openRegion286Bam <- gsub("\\.bam", "_openRegions286\\.bam", sortedBAM286)

library(rtracklayer)
export(atacReads286_Open, openRegion286Bam, format = "bam")

openRegion286BigWig <- gsub("\\.bam", "_openRegions286\\.bw", sortedBAM286)
openRegion286RPMBigWig <- gsub("\\.bam", "_openRegions286RPM\\.bw", sortedBAM286)
atacFragments286_Open <- granges(atacReads286_Open)
export.bw(coverage(atacFragments286_Open), openRegion286BigWig)


#annotating peakes to genes
library(ChIPseeker)
peaks286 <- "Peak Calling/SRR13169286.broad_peaks.broadPeakdiff2"
peaks287 <- "Peak Calling/SRR13169287.broad_peaks.broadPeakdiff2"
peaks288 <- "Peak Calling/SRR13169288.broad_peaks.broadPeakdiff2"
peaks289 <- "Peak Calling/SRR13169289.broad_peaks.broadPeakdiff2"

MacsCalls286_Anno <- annotatePeak(peaks286, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
MacsCalls287_Anno <- annotatePeak(peaks287, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
MacsCalls288_Anno <- annotatePeak(peaks288, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
MacsCalls289_Anno <- annotatePeak(peaks289, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
MacsCalls286_Anno

plotAnnoPie(MacsCalls286_Anno)
plotAnnoPie(MacsCalls287_Anno)
plotAnnoPie(MacsCalls288_Anno)
plotAnnoPie(MacsCalls289_Anno)


#GREAT analysis
library(rGREAT)

great_job286 <- submitGreatJob(peaks286, species = "mm10")
availableCategories(great_job286)
great_ResultTable286 = getEnrichmentTables(great_job286, category = "GO")
names(great_ResultTable286)
great_ResultTable286[["GO Biological Process"]][1:4]
save(great_ResultTable286, file= "Great_Results286.RData")

great_job287 <- submitGreatJob(peaks287, species = "mm10")
availableCategories(great_job287)
great_ResultTable287 = getEnrichmentTables(great_job287, category = "GO")
names(great_ResultTable287)
great_ResultTable287[["GO Biological Process"]][1:4]
save(great_ResultTable287, file= "Great_Results287.RData")

great_job288 <- submitGreatJob(peaks288, species = "mm10")
availableCategories(great_job288)
great_ResultTable288 = getEnrichmentTables(great_job288, category = "GO")
names(great_ResultTable288)
great_ResultTable288[["GO Biological Process"]][1:4]
save(great_ResultTable288, file= "Great_Results288.RData")

great_job289 <- submitGreatJob(peaks289, species = "mm10")
availableCategories(great_job289)
great_ResultTable289 = getEnrichmentTables(great_job289, category = "GO")
names(great_ResultTable289)
great_ResultTable289[["GO Biological Process"]][1:4]
save(great_ResultTable289, file= "Great_Results289.RData")
x = great_ResultTable289[["GO Biological Process"]][1:4]

great_table <- read.table("greatExportAll.tsv", sep = "\t")
great_table286 <- read.table("greatExportAll286.tsv", sep = "\t")
par(mar=c(4,15,1,1))

wrap.labels <- function(x,len){
  sapply(x, function(y) paste(strwrap(y,len), collapse = "\n"), USE.NAMES = FALSE)
}

barplot(great_table$V6[1:5], horiz = T, 
        xlab = 'FDR q value',
        names.arg = wrap.labels(great_table$V3[1:5],6),las =2,font.axis = 1,
        main = "GREAT GO Biological Process for SRR13169289", cex.names = 0.5)

barplot(great_table286$V6[1:5], horiz = T, 
        xlab = 'FDR q value',
        names.arg = wrap.labels(great_table286$V3[1:5],6),las =2,font.axis = 1,
        main = "GREAT GO Biological Process for SRR13169286", cex.names = 0.5)

#diff ATAC seq analysis
library(ChIPQC)

peaks <- dir(pattern = ".broadPeakdiff2", full.names = TRUE)
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple=TRUE)
names(myPeaks) = c("psDC1", "psDC2", "nsDC1", "nsDC2")
Group = factor(c("psDC", "psDC", "nsDC", "nsDC"))
consensusToCount = soGGi::runConsensusRegions(GRangesList(myPeaks), "none")

blklist = import.bed("file.bed.gz")

library(DiffBind)

library(Rsubread)

Counts286 = featureCounts(sortedBAM286, isPairedEnd = TRUE)

occurances = elementMetadata(Counts286) %>% as.data.frame %>% dplyr::select(-conensusIDs) %>%
  rowSums

table(occurances) %>% rev %>% cumsum


# Differential Peak Calling -----------------------------------------------
setwd("/class/projects/bioinf545-2021/group-1/ATACseq")

library(ggplot2)
library(reshape2) 
library(pander)
library(Hmisc)
library(pastecs)
library(DESeq2)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(pheatmap)

# org.Mm.eg.db broke this week, workaround until fixed:
# install(remotes)
#remotes::install_version("RSQLite", version = "2.2.5")
# BiocManager::install('org.Mm.eg.db')
# restart R 
#library(org.Mm.eg.db)

library(Rsubread)
#counts =featureCounts("mergedBam.bam", isPairedEnd = TRUE)
#counts = counts$counts
x = featureCounts("SRR13169286.pruned.bam", isPairedEnd = TRUE)
info286 = featureCounts("SRR13169286.pruned.bam", isPairedEnd = TRUE)$annotation    
info287 = featureCounts("SRR13169287.pruned.bam", isPairedEnd = TRUE)$annotation 
info288 = featureCounts("SRR13169288.pruned.bam", isPairedEnd = TRUE)$annotation 
info289 = featureCounts("SRR13169289.pruned.bam", isPairedEnd = TRUE)$annotation 
counts286 = featureCounts("SRR13169286.pruned.bam", isPairedEnd = TRUE)$counts
counts287 = featureCounts("SRR13169287.pruned.bam", isPairedEnd = TRUE)$counts
counts288 = featureCounts("SRR13169288.pruned.bam", isPairedEnd = TRUE)$counts
counts289 = featureCounts("SRR13169289.pruned.bam", isPairedEnd = TRUE)$counts

counts = cbind(counts286, counts287, counts288, counts289)
info <- read.table(file = "Peak Calling/SRR13169286.broad_treat_pileup.bdg", sep = "\t", header = TRUE)

df = as.data.frame(counts)
df = data.frame(info286$Chr, info286$Start, info286$End, counts286, counts287, counts286, counts289)
SampleID <- colnames(counts)
Treatment <- c("Stimulated", "Stimulated", "Not Stimulated", "Not Stimulated")
si = data.frame(SampleID, Treatment)
rownames(si) = si$SampleID
si
#colnames(si) = c("sample", "treatment")
#rownames(si) = si$sample 
#si$donor = factor(si$donor)

pander(dim(df), "Data dimensions") 
pander(head(df))
pander(head(si))

print(summary(si))

### Remove missing peaks 
df = df[apply(df[,1:ncol(df)], 1, max) > 50,]
pander(dim(df), "Data dimensions")

cm = df[,1:ncol(df)]
pander(quantile(rowSums(cm)))
pander(quantile(rowMeans(cm)))
pander(quantile(apply(cm,1,max)))

#### Data Exploration 
print(mean(df$SRR13169286.pruned.bam))

print(sd(df$SRR13169286.pruned.bam))

sim = rnorm(1000, mean=100, sd=10)
boxplot(sim)
hist(sim)

pander(quantile(sim))

# label outliers 
pw = df$SRR13169286.pruned.bam
uq = quantile(pw, 0.75)
print(mean(pw > 1.5 * uq))

## [1] 0.1635487

iqr = IQR(pw)
print(mean(pw > 3 * iqr))

## [1] 0.07711043

pw = log(df$SRR13169286.pruned.bam + 1)
boxplot(pw)
hist(pw)

# check for heteroscedasticity 
rowsummary = data.frame(rowmeans = apply(df[, 1:ncol(df)], 1, mean), rowsds = apply(df[, 1:ncol(df)], 1, sd))
ggplot(data=rowsummary, aes(x=rowmeans, y=rowsds)) + geom_point() + xlab("Peak means") + ylab("Peak SDs")

#### Data Normalization 
counts = df[,1:ncol(df)]
dds = DESeqDataSetFromMatrix(countData = counts[,order(colnames(counts))], colData = si ,design = ~ Treatment)
dds = DESeq(dds)

cm = data.frame(counts(dds, normalized=TRUE))
rownames(cm) = paste0(info$chr, '_', info$Start, '_', info$End)

#### Data Visualization 
lf = melt(cm, id.vars=c())
pander(head(lf))

ggplot(data=lf, aes(x=variable, y=value)) + geom_boxplot(aes(group=variable)) + xlab("Sample") + ylab("Normalized Count") + coord_flip()

ggplot(data=lf, aes(x=value)) + geom_freqpoly(aes(group=variable, color=variable), bins=30) + xlab("Sample") + ylab("Normalized Count")

libsize = data.frame(x=sizeFactors(dds), y=colSums(assay(dds)))
ggplot(data=libsize, aes(x=x, y=y)) + geom_point() + geom_smooth(method="lm") + xlab("Estimated size factor") + ylab("Library size")

#### PCA
pca = prcomp(t(cm))
print(summary(pca))
pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
pcaData=merge(pcaData, si)
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))
p=ggplot(data=pcaData, aes(x = PC1, y = PC2, color=celltype)) + geom_point(size=3)
p=p+xlab(paste0("PC1: ", percentVar[1], "% variance"))
p=p+ylab(paste0("PC2: ", percentVar[2], "% variance"))
print(p)

q=ggplot(data=pcaData, aes(x = PC3, y = PC4, color=celltype)) + geom_point(size=3)
q=q+xlab(paste0("PC3: ", percentVar[3], "% variance"))
q=q+ylab(paste0("PC4: ", percentVar[4], "% variance"))
print(q)

varexp = data.frame(x=1:length(percentVar), y=percentVar)
varexp$x = factor(varexp$x)
ggplot(data=varexp, aes(x=x, y=y)) + geom_bar(stat="identity") + xlab("Principal Component") + ylab("Proportion of variation (%)")

loadings = abs(pca$rotation)
contribution = as.data.frame(sweep(loadings, 2, colSums(loadings), "/"))
contribution = contribution[with(contribution, order(-PC1)),]
pander(head(contribution))

#### Annotate Genome Context 
gr = makeGRangesFromDataFrame(df, keep.extra.columns=T)
peakAnno = annotatePeak(gr, tssRegion=c(-1000, 1000), TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)

dfPA = as.data.frame(peakAnno)
rownames(dfPA) = paste0(dfPA$seqnames, '_', dfPA$start, '_', dfPA$end)
selpeaks = dfPA[rownames(head(contribution, 500)),]
pathway1 = enrichPathway(selpeaks[abs(selpeaks$distance) < 5000,]$geneId)
pander(head(pathway1))

dotplot(pathway1)

#### Differential Peak Calling
res = results(dds, lfcThreshold=1, contrast=c("Treatment", "Stimulated", "Not Treatment"))
print(mcols(res, use.names=T))

hist(res$pvalue, breaks=0:20/20, col="grey50", border="white", xlim=c(0,1), main="Histogram of p-values", xlab="p-value")
plotMA(res, ylim = c(-5, 5)) # visualize log fold changes 

print(sum(res$padj < 0.01 & abs(res$log2FoldChange) > 1))

mat = cm[which(res$padj < 0.01 & abs(res$log2FoldChange) > 1),]
mat = mat - rowMeans(mat)
anno = as.data.frame(colData(dds)[, c("sample", "celltype")])
rownames(mat) = NULL
pheatmap(mat, annotation_col = anno, scale="row")


selpeaks = dfPA[rownames(cm[which(res$padj < 0.1 & res$log2FoldChange>0),]),]
pathwayUp = enrichPathway(selpeaks[abs(selpeaks$distance) < 5000,]$geneId)
pander(head(pathwayUp))

selpeaks = dfPA[rownames(cm[which(res$padj < 0.1 & res$log2FoldChange<0),]),]
pathwayDown = enrichPathway(selpeaks[abs(selpeaks$distance) < 5000,]$geneId)
pander(head(pathwayDown))





