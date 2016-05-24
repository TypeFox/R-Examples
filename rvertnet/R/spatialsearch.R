#' Find records within some distance of a point given latitude and longitude.
#' 
#' Searches by decimal latitude and longitude to return any occurrence record
#' within the input distance (radius) of the input point.
#'
#' @export
#' @param lat Latitude of the central point, in decimal degrees (numeric) Required.
#' @param long Longitude of the central point, in decimal degrees (numeric) Required.
#' @param radius Radius to search, in meters (numeric). There is no default value for
#' this parameter. Required.
#' @param limit (numeric) Limit on the number of records returned. If >1000 results, we use
#' a cursor internally, but you should still get up to the results you asked for. See also 
#' \code{\link{bigsearch}} to get larger result sets in a text file via email.
#' @param compact Return a compact data frame (logical)
#' @param verbose Print progress and information messages. Default: TRUE
#' @param ... Curl arguments passed on to \code{\link[httr]{GET}}
#' @details \code{\link{spatialsearch}} finds all records of any taxa having decimal lat/long
#'    coordinates within a given radius (in meters) of your coordinates.
#' @return A data frame of search results
#' @references \url{https://github.com/VertNet/webapp/wiki/The-API-search-function}
#' @examples \dontrun{
#' res <- spatialsearch(lat = 33.529, long = -105.694, radius = 2000, limit = 10)
#' 
#' # Pass in curl options for curl debugging
#' library("httr")
#' out <- spatialsearch(lat = 33.529, long = -105.694, radius = 2000, limit = 10, config=verbose())
#' }

spatialsearch <- function(lat, long, radius, limit = 1000, compact = TRUE, verbose = TRUE, ...){
  args <- list(lat = lat, long = long, radius = radius)
  vertwrapper(fxn = "spatialsearch", args = args, lim = limit, compact = compact, verbose = verbose, ...)
}
