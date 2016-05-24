GetCommonNames <- function(MyEOLs, output=c("detail", "counts")) {
  MyEOLs <- RemoveNAFiles(MyEOLs)
  output <- match.arg(output)
  CommonNames <- matrix(nrow=0, ncol=4)
  colnames(CommonNames) <- c("Taxon", "eolID", "Common Name", "language")
  CNOverview <- data.frame(matrix(nrow=length(MyEOLs), ncol=2))
  colnames(CNOverview) <- c("Taxon", "eolID")
  for(i in sequence(length(MyEOLs))) {
    taxon <- NA
    eolID <- NA
    res <- PageProcessing(MyEOLs[i])$taxonConcept
    if(!is.null(res)) {
      taxon  <- res[[which(names(res) == grep("ScientificName", names(res), ignore.case=TRUE, value=T))]] #because some are cap and some are not
      eolID <- res$taxonConceptID
      CNOverview[i,1:2] <- c(taxon, eolID)
      CNs <- which(names(res) == "commonName")
      taxonCommonNames <- rep(NA, 4)
      for(j in sequence(length(CNs))) {
        if(any(names(res[[CNs[j]]]) == ".attrs")){
          language <- as.character(res[[CNs[j]]]$.attr[which(names(res[[CNs[j]]]$.attr) == "lang")]) #not tidy, but effective for multiple entries
          taxonCommonNames <- c(taxon, eolID, res[[CNs[j]]]$text, language)
          CommonNames <- rbind(CommonNames, taxonCommonNames, deparse.level=0)
          CommonNames <- data.frame(CommonNames, stringsAsFactors=FALSE)
          if(sum(grepl(language, colnames(CNOverview))) == 0) {
            CNOverview <- cbind(CNOverview, rep(0, length(MyEOLs)))
            colnames(CNOverview) <- append(colnames(CNOverview[-dim(CNOverview)[2]]), language)
          }
        }
        languageColumn <- which(colnames(CNOverview) == language) 
        CNOverview[i,languageColumn] <- as.numeric(CNOverview[i,languageColumn])+1   
      }
    }
  }
  if (output == "detail")
    return(CommonNames)
  if (output == "counts")
    return(CNOverview)
}