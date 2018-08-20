####Script for Leila for DE
source("http://www.Bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite( c("ShortRead","DESeq2", "edgeR") )
install.packages("statmod")
install.packages("corrplot")
library("edgeR")
library("corrplot")

#change this to correct path
setwd("/Users/srs57/Desktop/Jim/R")

#read in annotations
annot <- read.delim("coffea_canephora_functions2.csv", row.names=1, sep = "\t")

#change this to correct path
setwd("/Users/srs57/Desktop/Jim/R")

#Get counts with featureCount subread package, use as input
##Read in files, group, and format
y <- read.delim("gene_count_matrix13.csv", row.names=1, stringsAsFactors=FALSE, sep = ",")
raw_counts = y

#Filter
keep <- rowSums(cpm(y) > 5) >= 3  #keep genes with at least 5 cpm in at least 3  samples
y <- y[keep,]
table(keep)
dim(keep)

#Get fpkms
expressioncount <- read.delim("cds_len.csv", row.names=1)
result <- rpkm(y, normalized.lib.sizes=FALSE, gene.length=expressioncount$length)
write.csv(result, file="rpkm.csv")

#make object
samples <- read.delim("treatment.csv", row.names="Sample", header=TRUE, sep = ",")
samples <- cbind(samples, (make.names(paste(samples$Treatment, samples$Time, sep = ""))))
colnames(samples) <- c("Treatment", "Time", "Group")
group <- factor(samples$Group)
y <- DGEList(counts=y, group=samples$group)

##Re-compute the library sizes:
y$samples$lib.size <- colSums(y$counts)

#check library sizes
png(filename = 'counts.png', width = 1500, height = 1500, units = 'px')
barplot(y$samples$lib.size, main="Raw Counts", xlab = "sample", ylab = "counts")
dev.off()

# Get log2 counts per million
logcpm <- cpm(y$counts,log=TRUE, prior.count = 1)

# Check distributions of samples using boxplots
# Let's add a blue horizontal line that corresponds to the median logCPM
#color by group
par(mfrow=c(1,1),oma=c(2,0,0,0))
group.col <- c("red","blue", "green", "purple", "yellow", "orange", "skyblue","red","blue", "green", "purple", "yellow", "orange", "skyblue", "red", "blue", "green", "purple")[group]
png(filename = 'logcpm.png', width = 1500, height = 1500, units = 'px')
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,col=group.col,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs\n(coloured by groups)",cex.main=0.8)
dev.off()

#Estimate normalization factors, By default, calcNormFactors uses the TMM method and the sample whose 75%-ile (of library-scale-scaled counts) is closest to the mean of 75%-iles as the reference.
y = calcNormFactors(y)
y$samples

logcounts <- cpm(y,log=TRUE)
par(mfrow=c(1,2))
plotMD(logcounts,column=1)
abline(h=0,col="grey")
plotMD(y,column = 2)
abline(h=0,col="grey")
par(mfrow=c(1,1))

#estimate dipsersion
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
plotBCV(y)
summary(y$prior.df)
sqrt(y$common.disp)  #The square root of the common dispersion gives the coefficient of variation of biological variation (BCV).
png(filename = 'bcv.png', width = 800, height = 800, units = 'px')
plotBCV(y)
dev.off()

########QC############
png(filename = 'mds.png', width = 800, height = 800, units = 'px')
plotMDS(y, method="bcv")
dev.off()

png(filename = 'dendrogram.png', width = 800, height = 800, units = 'px')
tcounts <- t(as.table(as.matrix(y)))
counts.dist = hclust(dist(tcounts))
plot(counts.dist)
dev.off()

#scatterplot of reps
#png(filename = 'pnp1_pnp2_scatter.png', width = 800, height = 800, units = 'px')
y.cpm <-cpm(y, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)
#plot(y.cpm[,"M82_5d_1"], y.cpm[,"M82_5d_2"], log="xy")
#dev.off()

write.csv(y.cpm, file="norm_cpm_all.csv")

#Correlation Matrix
png(filename = 'corr_matrix.png', width = 1500, height = 1500, units = 'px')
corrplot(cor(y.cpm), method="square", cl.lim=c(0,1), tl.col="black", addgrid.col="black", is.corr=FALSE, number.cex = 0.5)
dev.off()

#############DE#######################
#Set up design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- glmFit(y, design)
logFC <- predFC(y,design,prior.count=1,dispersion=0.05)
cor(logFC[,1:6])  #correlation matrix of the pooled samples

#Estimate parameters
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y,design)

#Set up comparisons
# Get all genes that are de in cfm to wt comparisons
my.contrasts <- makeContrasts(
  C0T1 = CtrlT0 - TestT1,     #these should be switched, as written these will me it is downreg or upreg in control
  C0T3 = CtrlT0 - TestT3,
  C0T5 = CtrlT0 - TestT5,
  levels=design
)

contrast.names <- colnames(my.contrasts)
d <- data.frame(matrix(, nrow=nrow(y$counts), ncol=0))

#####To find genes diff:
prefix.lrt  <- glmLRT(fit, contrast=my.contrasts)

#loop through contrasts
for (i in c(1:length(contrast.names))){
  
  #####To find genes diff:
  prefix.lrt  <- glmLRT(fit, contrast=my.contrasts[,contrast.names[i]])
  
  # total number of genes significantly up-regulated or down-regulated at 5% FDR
  summary(prefix.dt <- decideTestsDGE(prefix.lrt, adjust.method="fdr", p.value=0.05, lfc=2))
  
  #toptags
  result <- topTags(prefix.lrt, n=Inf)
  keep <- result$table$FDR <= 0.05 
  keep_up <- result$table$logFC >= 1
  keep_down <- result$table$logFC <= -1
  #result2 <- result[keep,]
  result_up <- result[keep_up,]
  result_down <- result[keep_down,]
  
  # plot all the logFCs against average count size, highlighting the DE genes
  outfile3 <- paste(contrast.names[i], "_fdr0.05logFC2.png", sep="")
  png(filename = outfile3, width = 1000, height = 1000, units = 'px')
       cols <- rep("black",nrow(d))
       cols[prefix.lrt$table$PValue<.05] <- "red"
       plotSmear(prefix.lrt, col=cols) #de.tags=result2$Row.names)
       abline(h=c(-1,1), col="blue") #The blue lines indicate 2-fold up or down, ie log2(2)=1
  dev.off()
  
  #merge functions file
  result <- merge(result, annot, by=0, all=FALSE)
  result_up <- merge(result_up, annot, by=0, all=FALSE)
  result_down <- merge(result_down, annot, by=0, all=FALSE)
  
  #export results
  outfile <- paste(contrast.names[i], "_all.csv", sep="")
  write.table(result, file=outfile, sep = "\t")  
  
  outfile2 <- paste(contrast.names[i], "_fdr0.05_up.csv", sep="")
  write.table(result_up, file=outfile2, sep = "\t")
  
  outfile4 <- paste(contrast.names[i], "_fdr0.05_down.csv", sep="")
  write.table(result_down, file=outfile4, sep = "\t")
}
