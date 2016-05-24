#' @title SotkanetIndicators
#' @description SotkanetIndicators retrieves Sotkanet data corresponding to a
#' specified data identifier from 
#' \url{http://www.sotkanet.fi/rest/1.1/indicators}
#'
#' @param id Dataset identifier
#' @param type output format ("table" or "raw")
#' @return sotkanet json query in selected output format
#'
#' @export
#' @references See citation("sotkanet") 
#' @author Einari Happonen. Maintainer: Louhos/Opasnet \email{louhos@@googlegroups.com}
#' @examples \dontrun{sotkanet.indicators <- SotkanetIndicators(type = "table")}
#' @keywords utilities
SotkanetIndicators <- function(id = NULL, type = "table")
{

  base.url <- base_url()	
  url <- paste(base.url, 'indicators', sep = "")

  if (!is.null(id))
    url <- paste(url, id, sep='/')

  res <- sotkanet.json_query(url)

  if (type == "table") {
    res <- SotkanetCollect(res, "indicator")
  }

  res
}


