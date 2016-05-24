#' @title get PMIDs using Journal names and Keywords
#' @description retrieve PMIDs (each PMID is 8 digits long) from PubMed for Specific Journal, Keywords and date.
#' @param keys     keywords
#' @param journal journal name
#' @param dFrom   start year
#' @param dTo     end year
#' @param n       max number of retrieved articles
#' @seealso \code{\link{getAbstracts}}
#' @seealso \code{\link{editPMIDs}}
#' @seealso \code{\link{getPMIDs}}
#' @export
#' @examples
#' # getPMIDsByKeyWords(keys="breast cancer", journal="science",dTo=2013)
#' 
#' # getPMIDsByKeyWords(keys="breast cancer", journal="science")
#' 
#' # getPMIDsByKeyWords(keys="breast cancer",dFrom=2012,dTo=2013)
#' 
#' # getPMIDsByKeyWords(journal="science",dFrom=2012,dTo=2013)
getPMIDsByKeyWords <-function(keys=NULL,journal=NULL,dFrom=NULL, dTo=NULL, n=10000){
  eSearch <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=" 
  
  if(is.null(keys) && is.null(journal)){
    stop("at most one of 'keys' and 'journal' can be NULL!")
  }
  
  kQ <- ""
  if(!is.null(keys)){
    kQ <- str_replace_all(keys, " ", "+")
  }
  
  jQ <- ""
  if(!is.null(journal)){
    jt <- str_replace_all(journal, " ", "+")
    jQ <- paste(jt, "[journal]", sep = "")
  }
  
  dQ <- ""
  if(!is.null(dFrom) && !is.null(dTo)){
    d1 <- paste(dFrom, dTo, sep = ":")
    dQ <- paste(d1, "[pdat]", sep = "")
  }else if(!is.null(dFrom) && is.null(dTo)){
    dQ <- paste(dFrom, "[pdat]", sep = "")
  }else if(is.null(dFrom) && !is.null(dTo)){
    dQ <- paste(dTo, "[pdat]", sep = "")
  }
  
  hlpQ1 <- kQ
  hlpQ2 <- hlpQ1
  hlpQ3 <- hlpQ2
  
  if(str_length(jQ)>0){
    hlpQ2 <- paste(hlpQ1, "AND", jQ, sep = "+")
    hlpQ3 <- hlpQ2
  }
  
  if (str_length(dQ) > 0){
    hlpQ3 <- paste(hlpQ2, "AND", dQ, sep = "+")
  } 
  
  rmQ <- paste("&retmax=", n, sep="")
  hlpQ4 <- paste(hlpQ3, rmQ, sep="")
  
  searchUrl <- paste(eSearch, hlpQ4, sep = "" )
  
  hlpURL <- getURL(searchUrl)
  
  doc <- xmlTreeParse(hlpURL, asText = TRUE)     
  IdlistHlp = xmlValue(doc[["doc"]][["eSearchResult"]][["IdList"]])
  
  if (length(IdlistHlp) > 0){
    Idlist <- substring(IdlistHlp, seq(1, nchar(IdlistHlp)-1, 8), seq(8, nchar(IdlistHlp), 8))
  }
  return(Idlist)
}
