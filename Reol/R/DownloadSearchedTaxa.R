APItaxon <- function(taxon) {
  taxon <- gsub("_", "+", taxon)
  taxon <- gsub(" ", "+", taxon)  #doesn't work
  return(taxon)
}


MatchTaxatoEOLID <- function(ListOfTaxa, exact=TRUE, ...){ 
  #Match a search taxon to an EOLID for downloading or storing
  #API can support fuzzy matching up to 30 matches if exact=F, but then pages have to be specified.  Might be a good thing to add later. 
  eolPageNumbers <- rep(NA, length(ListOfTaxa))
  speciesNameForRef <- rep(NA, length(ListOfTaxa))
  for (i in sequence(length(ListOfTaxa))) {  
    taxon <- APItaxon(ListOfTaxa[i])
	web <- paste('http://eol.org/api/search/1.0.xml?q=', taxon, '&exact=', exact, '&page=1', sep="")
    a <- getURL(web, ...)
    searchRes <- NULL
    searchRes <- xmlToList(xmlRoot(xmlParse(a, getDTD=FALSE), ...), simplify=FALSE)
    if(searchRes$totalResults == 1) {  #didn't match any eol taxa
      eolPageNumbers[i] <- searchRes$entry$id  #there are other matches sometimes as well
      speciesNameForRef[i] <- searchRes$entry$title
    }
  }
  return(data.frame(ListOfTaxa, speciesNameForRef, eolPageNumbers, stringsAsFactors=F))
}


DownloadSearchedTaxa <- function(ListOfTaxa, to.file=TRUE, MyKey=NULL, exact=TRUE, verbose=TRUE, ...) {
  matches <- MatchTaxatoEOLID(ListOfTaxa, exact=exact, ...)
  pagesToDownload <- unique(matches[,3])
  fileNames <- data.frame(paste("eol", matches[,3], ".xml", sep=""), stringsAsFactors=FALSE)
  colnames(fileNames) <- "fileNames"
  if(any(!is.na(pagesToDownload))) {
    DownloadEOLpages(as.numeric(pagesToDownload), to.file, MyKey, verbose=verbose, ...)  
  }
}
