enrichment_test <-
function(myInterestingGenes, geneNames, geneID2GO, topNum) {
  # "myInterestingGenes" is the given set of genes
  # "geneNames" is the population of genes
  # "geneID2GO" is the mapping of gene IDs to GO terms
  #geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes),levels=c(0,1))
  names(geneList) <- geneNames
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata, classic = resultFis, orderBy = "weight", ranksOf = "classic", topNodes = topNum)
  return(allRes)
}
