PingAPI <- function(MyKey=NULL) {
  web <- "http://eol.org/api/ping.xml"
  if(!is.null(MyKey))
    web <- paste(web, "?key=", MyKey, sep="")
  a <- getURL(web)  
  xmlToList(xmlRoot(xmlParse(a, getDTD=FALSE)))$message
}