  #'Search biorxiv.org
#' @description return a list of URLS, IDS and number of results found from search of biorxiv.org
#' @param query the terms to search for in biorxiv
#' @param limit the number of results to return
#' @examples \dontrun{
#'   ### Get search results
#'   bxEco <- bx_search("ecology",limit=20) 
#'   summary(bxEco)
#' }
#' @details This uses the generic search interface, therefore be aware that you'll have to do post download filtering if you want terms in a specific field
#' @return a list with the following elements: a vector of URL's for bioRxiv papers from the search terms,and the maximum number of results
#' @export
#' @import XML
#' @importFrom utils URLencode

bx_search <- function(query, limit = 10){
  base <- "http://www.biorxiv.org/search/"
  page <- "?page="
  URL <- paste(base,URLencode(query),page,"0",sep="") 
  pgRes <- htmlParse(URL)
  
  ## Get total number of results possible
  ## Be careful, this may be a fragile way to extract the information, because it just relies on H1 from the header
  maxRes <- xpathApply(pgRes,"//h1",xmlValue)
  ### Remove Commas
  maxRes <- gsub(",", "", maxRes, fixed = TRUE) 
  maxRes <- as.numeric(regmatches(maxRes[[1]],regexpr('?[0-9]+' ,maxRes[[1]])))
  
  ## Handle a situation with no results
  if(identical(maxRes, numeric(0))){
    return(list(data=NULL,found=0))
  }
  ## Get the total number of pages
  pgCtLst  <-  xpathApply(pgRes, "//li/a")
  lastIndex <- which(unlist(lapply(pgCtLst,function(x){grepl("Last",xmlValue(x))}))==TRUE)
  ## Handle if there is just 1 page of results
  if(identical(lastIndex, integer(0))){
    lastPg <- 1
  } else {
    lastPg <- as.numeric(strsplit(xmlGetAttr(pgCtLst[[lastIndex]],"href"),"=")[[1]][2])
  }
  ## Get first round of URL's 
  ftURL <- unlist(lapply(xpathApply(pgRes, "//a[@class]",xmlGetAttr, "href"),grep,pattern="/content/early",value=T))
  ftURL <- sapply(ftURL,function(x){return(paste("http://www.biorxiv.org",x,sep=""))})
  
  ### Adjust for limits
  
  if(limit < maxRes){
    lastPg <- (limit%%10)+1
  }
  
  ### Loop through all the pages if more than one.
  
  if(lastPg > 0){
    for(i in 1:lastPg){
      URL <- paste(base,URLencode(query),page,i,sep="") 
      pgRes <- htmlParse(URL)
      
      tmpURL <-  unlist(lapply(xpathApply(pgRes, "//a[@class]",xmlGetAttr, "href"),grep,pattern="/content/early",value=T))
      tmpURL <- sapply(tmpURL,function(x){return(paste("http://www.biorxiv.org",x,sep=""))})
      ftURL <- c(ftURL,tmpURL)
    }
  } 
  if(limit < length(ftURL)){
    ftURL <- ftURL[1:limit]
  }
  ftURL <- as.matrix(ftURL)
  colnames(ftURL) <- "URL"
  IDs <- sapply(ftURL, function(m) {m1 <- strsplit(m,"/")[[1]]
                paste(m1[which(m1 =="early"):length(m1)],collapse="/")}
                )
  
  bioX <- structure(list(URL = ftURL,ID = unname(IDs), found = maxRes,query=query,limit=limit), class = "biorxiv_search")
  return(bioX)
  
}
#' Summary of search results
#' @description create a summary of search results
#' @param object the search object to create a summary of
#' @param ... extra parameters to pass
#' @export
#' 
summary.biorxiv_search <- function(object,...){
  class(object) <- "summary.bxso"
  return(object)
}

#'Print summary results
#'@description print a summary of search results
#'@param x the biorxiv search object to print
#'@param ... extra parameters to the print function
#'@export

`print.summary.biorxiv_search` <- function(x,...){
cat("Search term:", x$query,"\n")
cat("Number of results returned:",x$limit,"\n")
cat("Number of results found:",x$found,"\n")
}



