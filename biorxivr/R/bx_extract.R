#' Extract data from search results
#' @description Extract information about each paper in search results, returning author info, abstract and download metrics
#' @param bxso search results from bx_search()
#' @examples \dontrun{
#' bxEco <- bx_search("ecology",limit=5) 
#' bxEcoData <- bx_extract(bxEco)
#' ## See a plot of abstract views
#' plot(bxEcoData[[1]],type="abs")
#' 
#' ## See a plot of downloads
#' plot(bxEcoData[[1]],type="dl)
#' 
#' }
#' @details The metrics field will return a month by month data frame with abstract views and pdf downloads
#' @return a list holding the details about all the papers from search results. Each element in the list is an S3 object that holds all the salient details about the paper
#' @export
#' @import XML
bx_extract <- function(bxso){
  
  p <- sapply(bxso$URL,bx_extract_single,simplify=F)
  return(p)
}


#'Extract data from a single record
#'@description Generate an S3 object that represents a single paper
#'@param bxso_url the URL of a single biorxiv paper
#'@return a single S3 object with all the salient details about a paper
#'@export

bx_extract_single <- function(bxso_url){
  pgRes <- htmlParse(bxso_url)
  meta.names <- xpathApply(pgRes, "//meta[@name]",xmlGetAttr,"name")
  meta.data <- xpathApply(pgRes, "//meta[@name]",xmlGetAttr,"content")# DC.Contributor"
  names(meta.data) <- unlist(meta.names)
  nlist <- unlist(meta.names)
  paper <- list()
  authors <- list()
authors$names <- unname(unlist(meta.data[which(nlist == "DC.Contributor")]))
authors$email <- unname(unlist(meta.data[which(nlist == "citation_author_email")]))
  paper$title <- unname(unlist(meta.data[which(nlist == "DC.Title")]))
  paper$abstract <- unname(unlist(meta.data[which(nlist == "DC.Description")]))
  paper$date <- unname(unlist(meta.data[which(nlist == "DC.Date")]))
  paper$DOI <- unname(unlist(meta.data[which(nlist == "DC.Identifier")]))
  paper$fulltext_url <- unname(unlist(meta.data[which(nlist == "citation_pdf_url")]))

  metricRes <- htmlParse(paste(bxso_url,".article-metrics"))
  ## This returns a value for each element and we'll convert it to a matrix
  metrics <- unlist(xpathApply(metricRes, "//tbody/tr/td",xmlValue))
  ## If there are no metrics, at this point  the metrics dataframe will be null, hence this if statement
if(!is.null(metrics)){
  metrics <- data.frame(matrix(metrics,ncol=3,nrow=length(metrics)/3, byrow=T))
  metrics[,3] <- as.numeric(as.character(metrics[,3]))
  metrics[,2] <- as.numeric(as.character(metrics[,2]))
metrics[,1] <- do.call(c,lapply(strsplit(as.character(metrics[,1])," "),function(x) {as.Date(strptime(paste(paste(x,collapse="-"),"01",sep="-"),"%b-%Y-%d"))}))
  colnames(metrics) <- c("date","Abstract","PDF")
}
  
    
return(structure(list(authors = authors, paper = paper, metrics = metrics),class = "biorxiv_paper"))
}


#' plot metric details for a paper
#' @description plot a summary of the views a paper has had
#' @param x the paper to plot a summary of
#' @param type the data to plot, 'abs' for abstract views, 'dl' for PDF downloads
#' @param ... extra parameters to pass
#' @export
#' @importFrom graphics plot lines
plot.biorxiv_paper <- function(x,type="abs",...){
  bxp <- x

  if(is.null(bxp$metrics)){          
    stop("There are no metrics available for this manuscript")
  }
  if(type=="abs"){
    plot(bxp$metrics$date,bxp$metrics$Abstract,main = "Number of Abstract views",xlab="Date",ylab="Number of views")
    lines(bxp$metrics$date,bxp$metrics$Abstract)
  }
  if(type=="dl"){
    plot(bxp$metrics$date,bxp$metrics$PDF,main = "Number of PDF downloads",xlab="Date",ylab="Number of downloads")
    lines(bxp$metrics$date,bxp$metrics$PDF)
  }
}
