#' As featurecollection
#'
#' @export
#' @param x Input
#' @details Helper function to make a FeatureCollection list object for use
#' in vizualizing, e.g., with \code{leaflet}
#' @examples
#' str <- 'MULTIPOINT ((100.000 3.101), (101.000 2.100), (3.140 2.180),
#' (31.140 6.180), (31.140 78.180))'
#' # wkt2geojson(str, fmt = 2) %>%
#' #   as_featurecollection() %>%
#' #   lawn::view()
as_featurecollection <- function(x) {
  if (!is(x, "geojson")) stop("Must be of class geojson", call. = FALSE)
  x <- unclass(x)
  if (!"properties" %in% names(x)) x$properties <- list()
  list(type = "FeatureCollection", features = list(x))
}
