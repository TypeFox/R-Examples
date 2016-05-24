#' Star Name 
#'
#' This will return infomation about the requested star, 
#' @title Search for star infomation
#' @param star a string of an existing star
#' @return JSON object with infomation about the queried star
#' @keywords Name
#' @examples
#' \dontrun{
#' library(RCurl)
#' library(RJSONIO)
#' ###Return Infomation about the Sun
#' starname("Sun")
#' ###Return Infomation about the star HIP1 HD224700 Gli
#' starname("HIP1 HD224700 Gli")
#' ###Return Infomation about all of the stars in the system
#' starname("")
#' }
#' @export
starname <- function (star) {
  internetcheck <- url.exists("http://star-api.herokuapp.com", timeout = 10)
  if( internetcheck != TRUE)
    stop('Hacktheuniverse or your internet connection is down')
  urldata <- paste('http://star-api.herokuapp.com/api/v1/stars/', URLencode(star), sep = "")
  data <- getURL(urldata)
  dataFrame <- RJSONIO::fromJSON(data)
  return (dataFrame)
}
