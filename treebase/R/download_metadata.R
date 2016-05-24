# file: metadata.R
# author: Carl Boettiger <cboettig@gmail.com>
# date: 11 May 2011
# Description: Implementation of the TreeBASE API in R 
# With input from Ducan Temple Lang



#' Download the metadata on treebase using the OAI-MPH interface
#' @param query a date in format yyyy-mm-dd
#' @param by return all data "until" that date, 
#'   "from" that date to current, or "all"
#' @param curl if calling in series many times, call getCurlHandle() first and 
#'  then pass the return value in here. Avoids repeated handshakes with server.
#' @details query must be#'  download_metadata(2010-01-01, by="until")
#'  all isn't a real query type, but will return all trees regardless of date
#' @examples \dontrun{
#' Near <- search_treebase("Near", "author", max_trees=1)
#'  metadata(Near[[1]]$S.id)
#' ## or manualy give a sudy id
#' metadata("2377")
#' 
#' ### get all trees from a certain depostition date forwards ##
#' m <- download_metadata("2009-01-01", by="until")
#' ## extract any metadata, e.g. publication date:
#' dates <- sapply(m, function(x) as.numeric(x$date))
#' hist(dates, main="TreeBase growth", xlab="Year")
#' 
#' ### show authors with most tree submissions in that date range 
#' authors <- sapply(m, function(x){
#'    index <- grep( "creator", names(x))
#'      x[index] 
#' })
#' a <- as.factor(unlist(authors))
#' head(summary(a))
#' 
#' ## Show growth of TreeBASE 
#' all <- download_metadata("", by="all")
#' dates <- sapply(all, function(x) as.numeric(x$date))
#' hist(dates, main="TreeBase growth", xlab="Year")
#' 
#' ## make a barplot submission volume by journals
#' journals <- sapply(all, function(x) x$publisher)
#' J <- tail(sort(table(as.factor(unlist(journals)))),5)
#' b<- barplot(as.numeric(J))
#' text(b, names(J), srt=70, pos=4, xpd=T)
#' }
#' @export
download_metadata <- function(query="", by=c("all", "until", "from"),
                            curl=getCurlHandle()){
  by = match.arg(by)
  oai_url <- "http://treebase.org/treebase-web/top/oai?verb=" 
  list_record <- "ListRecords&metadataPrefix=oai_dc&"
  midnight <- "T00:00:00Z"
  query <- paste(oai_url, list_record, by, "=", query, midnight, sep="")
  metadata_from_oai(query, curl=curl)
}


#' Search the dryad metadata archive
#' @param study.id the dryad identifier
#' @param curl if calling in series many times, call getCurlHandle() first and 
#'  then pass the return value in here. Avoids repeated handshakes with server.
#' @return a list object containing the study metadata
#' @examples \dontrun{
#'   dryad_metadata("10255/dryad.12")
#' }
#' @export
dryad_metadata <- function(study.id, curl=getCurlHandle()){
  oai_url <- "http://www.datadryad.org/oai/request?verb=" 
  get_record <- "GetRecord&metadataPrefix=oai_dc&identifier=" 
  query <- paste(oai_url, get_record, "oai:datadryad.org:", study.id, sep="")
  metadata_from_oai(query, curl=curl)
}





#' Get the metadata associated with the study in which the phylogeny
#'  was published.
#' @param study.id The treebase study id (numbers only, specify in quotes)
#' @param curl if calling in series many times, call getCurlHandle()
#'  first and then pass the return value in here.  avoids repeated
#' handshakes with server. 
#' @details if the tree is imported with search_treebase, 
#' then this is in tree$S.id
#' @keywords utilities internal
# @examples \dontrun{
#   tree <- search_treebase("1234", "id.tree")
#   show_metadata(tree$S.id)
# }
show_metadata <- function(study.id, curl=getCurlHandle()){
  oai_url <- "http://treebase.org/treebase-web/top/oai?verb=" 
  get_record <- "GetRecord&metadataPrefix=oai_dc&identifier=" 
  query <- paste(oai_url, get_record, "TB:s", study.id, sep="")
  out <- metadata_from_oai(query, curl=curl)
  out[[1]]
}


#' return the study.id from the search results.  
#'
#' get_study_id is deprecated, and now can be performed more easily using
#' phylo_metadata and oai_metadata search functions.  
#' @param search_results the output of download_metadata, or a subset thereof
#' @return the study id
#' @details this function is commonly used to get trees corresponding
#'   to the metadata search.  
# @examples \dontrun{
# all <- download_metadata("", by="all")
# 
# nature <- sapply(all, function(x) length(grep("Nature", x$publisher))>0)
# science <- sapply(all, function(x) length(grep("^Science$", x$publisher))>0)
# s <- get_study_id( all[nature] )
# }
#' @keywords internal
get_study_id <- function(search_results){
  sapply(search_results, 
          function(x){
            id <- x$identifier
            id <- sub(".*TB2:S(\\d+)+", "\\1", id)
          })
}


#' return the trees in treebase that correspond to the search results
#' get_study is deprecated, and now can be performed more easily using
#' phylo_metadata and oai_metadata search functions.  
#' @param search_results the output of download_metadata, or a subset thereof
#' @param curl the handle to the curl web utility for repeated calls, see
#'  the getCurlHandle() function in RCurl package for details.  
#' @param ... additional arguments to pass to search_treebase
#' @return all corresponding phylogenies.  
#' @details this function is commonly used to get trees corresponding
#'   to the metadata search.  
# @examples \dontrun{
# all <- download_metadata("", by="all")
# nature <- sapply(all, function(x) length(grep("Nature", x$publisher))>0)
# science <- sapply(all, function(x) length(grep("^Science$", x$publisher))>0)
# s <- get_study( all[nature] )
# s <- get_study(all[science])
# }
#' @keywords internal
get_study <- function(search_results, curl=getCurlHandle(), ...){
  sapply(search_results, function(x) search_treebase(x, input="id.study", curl=curl, ...))
}



#' Internal function for OAI-MPH interface to the Dryad database
#' @param query a properly formed url query to dryad
#' @param curl if calling in series many times, call getCurlHandle() first and 
#'  then pass the return value in here. Avoids repeated handshakes with server.
#' @keywords internal
#' @seealso \code{\link{dryad_metadata}}
metadata_from_oai <- function(query, curl=curl){
  message(query)
  tt <- getURLContent(query, followlocation=TRUE, curl=curl)
  doc <- xmlParse(tt)
  dc = getNodeSet(doc, "//dc:dc", namespaces=c(dc="http://www.openarchives.org/OAI/2.0/oai_dc/"))
  lapply(dc, xmlToList)
}


