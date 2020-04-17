setwd("/Users/srs57/Desktop/Zhilong/rQTL/linked_region2/")
mapthis <- read.cross(format=c("csv"), dir="/Users/srs57/Desktop/Zhilong/rQTL/linked_region2/", "table3", estimate.map=FALSE)

#drop markers with no data
mapthis<- drop.nullmarkers(mapthis)
plotMissing(mapthis)

#plot missing
png(filename = 'missing.png', width = 800, height = 800, units = 'px')
plotMissing(mapthis)
dev.off()

png(filename = 'geno_by_ind.png', width = 800, height = 800, units = 'px')
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
dev.off()

png(filename = 'geno_by_markers.png', width = 800, height = 800, units = 'px')
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")
dev.off()

#look for seg dist and drop bad markers
gt <- geno.table(mapthis)
gt[gt$P.value < 0.05/totmar(mapthis),]
todrop <- rownames(gt[gt$P.value < 1e-7,])
mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)

#fill gaps
mapthis <- fill.geno(mapthis, method=c("no_dbl_XO"), error.prob=0.01, map.function=c("haldane"))

#plot missing after filling gaps
png(filename = 'missing.png', width = 800, height = 800, units = 'px')
plotMissing(mapthis)
dev.off()

png(filename = 'geno_by_ind.png', width = 800, height = 800, units = 'px')
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
dev.off()

png(filename = 'geno_by_markers.png', width = 800, height = 800, units = 'px')
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")
dev.off()

#check no duplicate individuals
cg <- comparegeno(mapthis)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

#Get genotype frequencies
png(filename = 'geno_freq.png', width = 800, height = 800, units = 'px')
par(mfrow=c(1,3), las=1)
g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
for(i in 1:3) plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))
dev.off()
par(mfrow=c(1,1), las=1)

#Look at genotypes per chromosome
filenames <- "geno_per_chromo"
for(i in 1:12){
  outfile = paste(filenames, i, "png", sep=".")
  png(filename = outfile, width = 800, height = 800, units = 'px')
  plotGeno(mapthis, chr=i, ind=c(1:58))
  dev.off()
}

#estimate rf, takes long
mapthis <- est.rf(mapthis)   

#See if markers are switched
checkAlleles(mapthis)

#estimate map, takes awhile
newmap <- est.map(mapthis, error.prob=0.01, verbose=TRUE)
par(mfrow=c(1,1), las=1)

#plot rf
filenames <- "rf"
for(i in 1:12){
  outfile = paste(filenames, i, "png", sep=".")
  png(filename = outfile, width = 800, height = 800, units = 'px')
  plot.rf(mapthis, chr=i)
  dev.off()
}

#compare old map with new map
plotMap(mapthis, newmap)

#replace old map with new
hyper <- replace.map(mapthis, newmap)
par(mfrow=c(1,1), las=1)
plotMap(hyper, show.marker.names=TRUE, alternate.chrid=TRUE)  

#Find error LOD scores and remove ones above a cutoff, runs a long time
hyper <- calc.errorlod(hyper, error.prob=0.01)
print(toperr <- top.errorlod(hyper, cutoff=5))
bad <- toperr$marker
hyper<- drop.markers(hyper, bad)
summary(hyper)

#Plot crossovers
png(filename = 'crossovers.png', width = 800, height = 800, units = 'px')
nxo <- countXO(hyper)
dev.off()

png(filename = 'co_count.png', width = 800, height = 800, units = 'px')
plot(nxo, ylab="No crossovers")
dev.off()

#calculate QTL probabilities
hyper <- calc.genoprob(hyper, step=1, error.prob=0.1, stepwidth=c("fixed"), off.end=0)

#perform QTL genome scan with extended Haley-Knott regression
out.ehk <- scanone(hyper, method="hk")
png(filename = 'ehk.png', width = 800, height = 800, units = 'px')
plot(out.ehk, chr=c(1:12), ylim=c(0,5), gap=30, bandcol="gray90")
dev.off()

#IMpute
hyper <- sim.geno(hyper, step=1, n.draws=100, error.prob=0.1)
out2.imp <- scanone(hyper, method="imp")
tiff("whole_geno_imp.tif", res=200, compression = "lzw", height=3, width=6, units="in")
plot(out.imp, chr=c(1:12), ylim=c(0,5), gap=30, bandcol="gray90", xlab="Distance (cM)", ylab="LOD")
dev.off()

tiff("imp.tif", res=200, compression = "lzw", height=2.5, width=5, units="in")
plot(out.imp, chr=c(2), ylim=c(0,5), gap=30, bandcol="gray90", xlim=c(54950,55200),  xlab="Chromosome 2 (cM)", xaxt='n',ylab="LOD", pch=19, cex.lab=0.5)
dev.off()

#plot difference between imputed LOD and ehk LOD
png(filename = 'imp2ehk.png', width = 800, height = 800, units = 'px')
plot(out.imp - out2.em, chr=c(1:12), ylim=c(-0.5, 0.5), ylab=expression(LOD[IMP] - LOD[EM]))
abline(h=0, lty=3)
dev.off()

#find range of qtl
lodint(out.ehk, chr=2,  expandtomarkers=TRUE)
lodint(out.ehk, chr=8,  expandtomarkers=TRUE)

###Two QTL scan, http://www.rqtl.org/rqtltour2.pdf
hyper <- calc.genoprob(hyper, step=2)
out2 <- scantwo(hyper, method="hk")

png(filename = '2q_scan.png', width = 800, height = 800, units = 'px')
  plot(out2)
dev.off()

png(filename = '2q_scan2.png', width = 800, height = 800, units = 'px')
  plot(out2, lower="fv1")
dev.off()

png(filename = '2q_scan3.png', width = 800, height = 800, units = 'px')
  plot(out2, lower="fv1", upper="av1")
dev.off()

#permutation test, 2qtl
operm2 <- scantwo(hyper, method="hk", n.perm=5)
summary(out2, perms=operm2, alpha=0.2, pvalues=TRUE)

###find heritability of 2 QTL
#make QTL object
chrom <- c("2", "8") #make chromosome vector
coord <- c(55079.19, 50982.23) 
f2 <- subset(hyper, chr=chrom)
qtl <- makeqtl(hyper, chrom, coord, what="draws")

# fit model with 2 interacting QTLs interacting (performing a drop-one-term analysis)
#lod <- fitqtl(hyper, pheno.col=1, qtl, formula=y~Q1*Q2, method="hk")
#summary(lod)

# fit an additive QTL model
#lod.add <- fitqtl(hyper, pheno.col=1, qtl, formula=y~Q1+Q2, method="hk")
#summary(lod.add)

#write cross data
write.cross(hyper, format=c("csv"),filestem="PtoS_LA2109_linked_qtl.csv")

