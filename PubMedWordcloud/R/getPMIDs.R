#' @title get PMIDs using author names
#' @description retrieve PMIDs (each PMID is 8 digits long) from PubMed for author and the specified date.
#' @param author author's name
#' @param dFrom  start year
#' @param dTo    end year
#' @param n      max number of retrieved articles
#' @seealso \code{\link{getAbstracts}}
#' @seealso \code{\link{editPMIDs}}
#' @export
#' @examples
#' # getPMIDs(author="Yan-Hui Fan",dFrom=2007,dTo=2013,n=10)
#' 
#' # getPMIDs(author="Yanhui Fan",dFrom=2007,dTo=2013,n=10)
getPMIDs <-function(author,dFrom, dTo, n=50){
  eSearch <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=" 
  aL <- str_replace_all(author, " ", "+")
  aQ <- paste(aL, "[author]", sep = "")
  
  dQ <- ""
  if((str_length(dFrom) > 0) & (str_length(dTo) > 0)){
    d1 <- paste(dFrom, dTo, sep = ":")
    dQ <- paste(d1, "[pdat]", sep = "")
  }else if((str_length(dFrom) > 0) & (str_length(dTo) == 0)){
    dQ <- paste(dFrom, "[pdat]", sep = "")
  }else if((str_length(dTo) > 0) & (str_length(dFrom) == 0)){
    dQ <- paste(dTo, "[pdat]", sep = "")
  }
  
  hlpQ1 <- aQ  
  if (str_length(dQ) > 0){
    hlpQ1 <- paste(aQ, dQ, sep = "+")
  } 
  
  rmQ <- paste("&retmax=", n, sep="")
  hlpQ2 <- paste(hlpQ1, rmQ, sep="")
  
  searchUrl <- paste(eSearch, hlpQ2, sep = "" )
  hlpURL <- getURL(searchUrl)
  
  doc <- xmlTreeParse(hlpURL, asText = TRUE)     
  IdlistHlp = xmlValue(doc[["doc"]][["eSearchResult"]][["IdList"]])
  
  if (length(IdlistHlp) > 0){
    Idlist <- substring(IdlistHlp, seq(1, nchar(IdlistHlp)-1, 8), seq(8, nchar(IdlistHlp), 8))
  }
  return(Idlist)
}
