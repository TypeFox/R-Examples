#' Local Group Galaxy Name
#'
#' This will return infomation about the requested local group galaxy 
#' @title Search for local group galaxy infomation
#' @param localgroup a string of an existing local group galaxy
#' @return JSON object with infomation about the queried local group galaxy
#' @keywords Name
#' @examples
#' \dontrun{
#' library(RCurl)
#' library(RJSONIO)
#' ###Return Infomation about the local group galaxy 
#' localgroupname("IC 10")
#' ###Return Infomation about the local group galaxy 
#' localgroupname("WLM")
#' ###Return Infomation about all local group galaxies in the system
#' localgroupname("")
#' }
#' @export
localgroupname <- function (localgroup) {
  internetcheck <- url.exists("http://star-api.herokuapp.com", timeout = 10)
  if( internetcheck != TRUE)
    stop('Hacktheuniverse or your internet connection is down')
  urldata <- paste('http://star-api.herokuapp.com/api/v1/local_groups/', URLencode(localgroup), sep = "")
  data <- getURL(urldata)
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}
