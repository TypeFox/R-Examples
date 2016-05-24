#' Get districts for a latitude/longitude.
#'
#' @import httr
#' @template cg
#' @param latitude latitude of coordinate
#' @param longitude longitude of coordinate
#' @return List of output fields.
#' @export
#' @examples \dontrun{
#' cg_getdistrictlatlong(latitude = 35.778788, longitude = -78.787805)
#' }
cg_getdistrictlatlong <- function(latitude = NULL, longitude = NULL,
  key=getOption("SunlightLabsKey", stop("need an API key for Sunlight Labs")),
  callopts = list())
{
  url = "http://services.sunlightlabs.com/api/districts.getDistrictFromLatLong.json"
  args <- suncompact(list(apikey = key, latitude = latitude, longitude = longitude))
  tt <- GET(url, query=args, callopts)
  stop_for_status(tt)
  out <- content(tt, as = "text")
  res <- fromJSON(out, simplifyVector = FALSE)$response
  class(res) <- "cg_getdistrictlatlong"
  return( res )
}
