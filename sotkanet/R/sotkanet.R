#' Description:
#' SotkanetRegions retrieves Sotkanet regions data from
#' \url{http://www.sotkanet.fi/rest/1.1/regions}
#'
#' Arguments:
#'   @param type Return format ("table" or "raw")
#'
#' Returns:
#'   @return sotkanet json query in selected format
#'
#' @export
#' @references
#' See citation("sotkanet") 
#' @author Einari Happonen / Opasnet / Louhos. Maintainer: Louhos/Opasnet \email{louhos@@googlegroups.com}
#' @examples # sotkanet.regions <- SotkanetRegions(type = "table")
#' @keywords utilities

SotkanetRegions <- function(type = "table")
{

  base.url <- base_url()	
  url <- paste(base.url, 'regions', sep = "")

  res <- sotkanet.json_query(url)

  if (type == "table") {
    res <- SotkanetCollect(res, "region")
  }

  res
}


