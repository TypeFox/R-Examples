GetIUCNStat <- function(MyEOLs) {
  MyEOLs <- RemoveNAFiles(MyEOLs)
  IUCN <- matrix(ncol=3, nrow=length(MyEOLs))
  colnames(IUCN) <- c("Taxon", "eolID", "IUCNstat")
  for(i in sequence(length(MyEOLs))){
    res <- PageProcessing(MyEOLs[i])
    IUCNstat <- NA
    if(length(grep("IUCNConservationStatus", res)) > 0)
      IUCNstat <- res[[which(res == grep("IUCNConservationStatus", res, value=TRUE))]]$description
 scientificName  <- res$taxonConcept[[which(names(res$taxonConcept) == grep("ScientificName", names(res$taxonConcept), ignore.case=TRUE, value=TRUE))]] #because some are cap and some are not
  IUCN[i,] <- c(scientificName, res$taxonConcept$taxonConceptID, IUCNstat)
  }  
  return(data.frame(IUCN, stringsAsFactors=FALSE))
}