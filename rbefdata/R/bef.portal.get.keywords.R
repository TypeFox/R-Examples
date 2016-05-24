#' Fetch keywords from a BEFdata portal.
#'
#' This function fetches keywords from a BEFdata portal.
#'
#' @param curl You can pass in a curl handle with additional options. By default a curl handle is
#'        used to improve the memory footprint.
#' @param \dots Arguments passed to \code{\link[RCurl]{getURLContent}}
#'
#' @return The function returns a data frame of keyword (id, name, cound).
#'
#' @examples \dontrun{
#'             keywords=bef.portal.get.keywords()
#'           }
#' @import RCurl
#' @import XML
#' @export bef.portal.get.keywords bef.get.keywords
#' @aliases bef.get.keywords

bef.portal.get.keywords <- bef.get.keywords <- function(curl = getCurlHandle(), ...) {
   raw_keywords_xml = getURLContent(paste0(bef.options('url'),"/keywords.xml"), curl = curl, ...)
   if(getCurlInfo(curl)$response.code != 200) stop("Server Error. Try again later!")
   keywords_xml = xmlTreeParse(raw_keywords_xml, useInternalNodes = T)
   data_frame_keywords = xmlToDataFrame(keywords_xml, colClasses=c('numeric', 'character', 'numeric'), stringsAsFactors = FALSE)
   return(data_frame_keywords)
}
