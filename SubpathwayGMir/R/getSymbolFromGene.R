getSymbolFromGene <-
function(geneList){
    if(!exists("k2ri")) k2ri<-initializeK2ri()
      geneList <- as.character(geneList)
      gene2symbol <- GetK2riData("gene2symbol")
      symbolList  <- unique(as.character(gene2symbol[gene2symbol[,1] %in% geneList,2]))
      return(symbolList)
}
