#source("http://bioconductor.org/biocLite.R")
#biocLite(c("hgu95av2.db", "GO.db"))
#biocLite("goProfiles")
#biocLite("topGO")
#library("GO.db")
#library("goProfiles")
#library("AnnotationDbi")
#library(topGO)

#http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html

#for treatments with 2 genotypes, significant genes that were found in both genotype comparisons were selected first before running GO
setwd("~/Desktop/Jim/R")
fileNames <- Sys.glob("C*.txt")

for (file in fileNames) {
  file <- sub("\"", "", file)
  geneID2GO <- readMappings(file = "go_assoc.tsv")  
  geneUniverse <- names(geneID2GO)
  genesOfInterest <- read.delim(file,header=TRUE)
  genesOfInterest = genesOfInterest[-1,]
  genesOfInterest <- as.vector(genesOfInterest) 
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  #read in annotations
  annot <- read.delim("coffea_canephora_functions4.csv", sep = ",", header=TRUE)
  
  #make object
  myGOdata <- new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  #The list of genes of interest can be accessed using the method sigGenes():
  sg <- sigGenes(myGOdata)
  str(sg)
  numSigGenes(myGOdata) 
  
  #enrichment test
  resultTopgo <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
  # see how many results we get where weight01 gives a P-value <= 0.05:
  mysummary <- summary(attributes(resultTopgo)$score <= 0.05)
  if(nrow(as.matrix(mysummary)) == 3){
    numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.05
    allRes <- GenTable(myGOdata, classicFisher = resultTopgo, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = numsignif)
    
    #get genes with sig go
    myterms <- allRes$GO.ID
    mygenes <- genesInTerm(myGOdata, myterms)
    m <- data.frame(NULL)
    
    for (i in 1:length(myterms)){
      myterm <- myterms[i]
      mygenesforterm <- mygenes[myterm][[1]] 
      myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
      mygenesforterm2 <- mygenesforterm[myfactor == TRUE] 
      #print(paste("Term",myterm,"genes:",mygenesforterm2))
      if (length(mygenesforterm2 >0)){
        m <- rbind(m, data.frame(myterm, mygenesforterm2))
      }
    }  
    
    #write output
    m <- merge(x=m, y=annot, by.x='mygenesforterm2', by.y='id')
    result <- paste(file, "_cc_funct.csv", sep="")
    write.csv(m, file=result)
    result2 <- paste(file, "_cc_go.csv", sep="")
    write.csv(allRes, file = result2)
  }
}
