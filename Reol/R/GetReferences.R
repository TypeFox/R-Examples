GetReferences <- function(MyEOLs, output=c("detail", "counts")) {
  MyEOLs <- RemoveNAFiles(MyEOLs)
  output <- match.arg(output)
  ReferenceList <- matrix(nrow=0, ncol=3)
  colnames(ReferenceList) <- c("Taxon", "eolID", "Reference")
  RefCounts <- matrix(nrow=length(MyEOLs), ncol=3)
  colnames(RefCounts) <- c("Taxon", "eolID", "Number Of References")
  for(i in sequence(length(MyEOLs))) {
    res <- PageProcessing(MyEOLs[i])$taxonConcept
    whichReferences <- which(names(res) == "reference")
    scientificName  <- res[[which(names(res) == grep("ScientificName", names(res), ignore.case=TRUE, value=T))]] #because some are cap and some are not
    RefCounts[i,] <- c(scientificName, res$taxonConceptID, length(whichReferences))
    for(j in sequence(length(whichReferences))) {
      ReferenceList <- rbind(ReferenceList, c(scientificName, res$taxonConceptID, as.character (res[whichReferences[j]])))
      #}
    }
  }
  RefCounts <- data.frame(RefCounts, stringsAsFactors=F)
  ReferenceList <- data.frame(ReferenceList, stringsAsFactors=F)
  if(output == "detail")
    return(ReferenceList)
  if(output == "counts")
    return(RefCounts)
}
