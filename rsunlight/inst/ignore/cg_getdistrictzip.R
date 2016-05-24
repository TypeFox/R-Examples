#' Get districts that overlap for a certain zip code.
#'
#' @import httr
#' @param zip zip code to search
#' @template cg
#' @return List of output fields.
#' @export
#' @examples \dontrun{
#' cg_getdistrictzip(zip = 27511)
#' }
cg_getdistrictzip <- function(zip = NULL,
  key=getOption("SunlightLabsKey", stop("need an API key for Sunlight Labs")),
  callopts = list())
{
  url = "http://services.sunlightlabs.com/api/districts.getDistrictsFromZip.json"
  args <- suncompact(list(apikey = key, zip = zip))
  tt <- GET(url, query=args, callopts)
  stop_for_status(tt)
  out <- content(tt, as = "text")
  res <- fromJSON(out, simplifyVector = FALSE)$response
  class(res) <- "cg_getdistrictzip"
  return( res )
}
