# This script takes the importedGTF and creates the coverage of each gene (meaning a data frame with 1.genename 2. Start 3. End)

# THIS IS A SIMPLE, NAIVE IMPLEMENTATION, WRITE THIS LATER MORE EFFICIENT!

getGeneLocation <- function(gtf){
  gtf$gene_id <- as.character(gtf$gene_id)
  gtf$gene_id <- gsub(" ","",gtf$gene_id)
  geneIds <- unique(gtf$gene_id)
  result <- data.frame(geneID="0", chr="0", start=0, end=0)
  for(i in 1:length(geneIds)){
    temp <- gtf[is.element(gtf$gene_id,geneIds[i]),c(1,4:5)]
    rangeTemp <- range(temp[,c(2:3)])
    enterThis <- data.frame(geneID=geneIds[i], chr=as.character(temp[1,1]), start=rangeTemp[1], end=rangeTemp[2])
    result <- rbind(result, enterThis)
  }
  result[-1,]
}