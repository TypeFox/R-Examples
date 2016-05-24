#' Get districts for a latitude/longitude or zips
#'
#' @export
#' @param latitude (numeric) latitude of coordinate
#' @param longitude (numeric) longitude of coordinate
#' @param zip (integer) A 5 digit zip code
#' @template districts
#' @template cg_query
#' @return List including data.frame and metadata about results, list, or httr
#' response object.
#' @details A zip code may intersect multiple Congressional districts, so it is not as precise as
#' using a latitude and longitude. In general, we recommend against using zip codes to look up
#' members of Congress. For one, it's imprecise: a zip code can intersect multiple congressional
#' districts. More importantly, zip codes are not shapes. They are lines (delivery routes), and
#' treating them as shapes leads to inaccuracies.
#' @examples \dontrun{
#' cg_districts(zip = 27511)
#' cg_districts(latitude = 35.778788, longitude = -78.787805)
#'
#' # most parameters are vectorized, pass in more than one value
#' cg_districts(zip = c(27511, 97202))
#' }

cg_districts <- function(latitude = NULL, longitude = NULL, zip = NULL, query=NULL,
  per_page=20, page=1, order = NULL, key = NULL, as = 'table', ...) {

  key <- check_key(key)
  args <- sc(list(apikey = key, latitude = latitude, longitude = longitude, zip = zip,
                          query = query, per_page = per_page, page = page, order = order))

  give_cg(as, cgurl(), "/districts/locate", args, ...)
}
