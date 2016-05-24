#' Download PDF's
#' @description Download PDF's of all the papers in your search results.
#' @param bxso search results from bx_search()
#' @param directory The location you want to download the PDF's to.
#' @param create TRUE or FALSE. If true create the directory if it doesn't exist
#' @examples \dontrun{
#' bxEco <- bx_search("ecology",limit=20) 
#' bx_download(bxEco,"~/biorxivPDF)
#' }
#' @export
#' @importFrom utils download.file


bx_download <- function(bxso, directory, create = TRUE){
  if(!file.exists(directory)  && create){
    dir.create(file.path(directory))
  }
  
  if(substr(directory,nchar(directory),nchar(directory)) != "/"){
    directory <- paste(directory,"/",sep="")
  }
  
  for(i in 1:length(bxso$URL)){
    URL <- paste(bxso$URL[i],".full.pdf",sep="") 
    download.file(URL,destfile = paste(directory,strsplit(bxso$ID[i],"/")[[1]][5],".pdf",sep="") ,method = "auto")
  }
  
}