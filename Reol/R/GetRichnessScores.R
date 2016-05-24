GetRichnessScores <- function(MyEOLs) {
  MyEOLs <- RemoveNAFiles(MyEOLs)
  richnessDF <- matrix(nrow=length(MyEOLs), ncol=3)
  for(i in sequence(length(MyEOLs))) {
    richnessData <- rep(NA, 3)
    res <- PageProcessing(MyEOLs[i])$taxonConcept
    scientificName  <- res[[which(names(res) == grep("ScientificName", names(res), ignore.case=TRUE, value=T))]] #because some are cap and some are not
    richnessData <- c(scientificName, res$taxonConceptID, res$additionalInformation$richness_score)
    richnessDF[i,] <- richnessData
  }
  richnessDF <- as.data.frame(richnessDF, stringsAsFactors=FALSE)
  colnames(richnessDF) <- c("Taxon", "eolID", "Richness_Score")
  return(richnessDF)
}
