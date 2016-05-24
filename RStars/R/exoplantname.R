#' Exoplanet Name
#'
#' This will return infomation about the requested exoplant, 
#' @title Search for exoplant infomation
#' @param exoplant a string of an existing exoplant
#' @return JSON object with infomation about the queried exoplant
#' @keywords Name
#' @examples
#' \dontrun{
#' library(RCurl)
#' library(RJSONIO)
#' ###Return Infomation about the exoplanet 11 Com
#' exoplantname("11 Com")
#' ###Return Infomation about the exoplanet Kepler-9
#' exoplantname("Kepler-9")
#' ###Return Infomation about all exoplants in the system
#' exoplantname("")
#' }
#' @export

exoplantname <- function (exoplant) {
  internetcheck <- url.exists("http://star-api.herokuapp.com", timeout = 10)
  if( internetcheck != TRUE)
    stop('Hacktheuniverse or your internet connection is down')
  urldata <- paste('http://star-api.herokuapp.com/api/v1/exo_planets/', URLencode(exoplant), sep = "")
  data <- getURL(urldata)
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}
