#' Visualize geojson from a character string or list
#'
#' @export
#' @param x Input, a geojson character string or list
#' @param center (numeric) A length two vector of the form: \code{longitude, latitude}
#' @param zoom (integer) A number between 1 and 18 (1 zoomed out, 18 zoomed in)
#' @param fmt Number of digits to display after the decimal point when formatting
#' coordinates.
#' @seealso \code{\link{as_featurecollection}}
#' @return Opens a map with the geojson object(s) using \code{leaflet}
#' @examples \dontrun{
#' # point
#' str <- "POINT (-116.4000000000000057 45.2000000000000028)"
#' wktview(str)
#'
#' # multipoint
#' df <- us_cities[1:5,c('long','lat')]
#' str <- multipoint(df)
#' wktview(str)
#' wktview(str, center = c(-100,40))
#' wktview(str, zoom = 3)
#'
#' # linestring
#' wktview(linestring(c(100.000, 0.000), c(101.000, 1.000), fmt=2))
#'
#' # polygon
#' a <- polygon(list(c(100.001, 0.001), c(101.12345, 0.001), c(101.001, 1.001), c(100.001, 0.001)))
#' wktview(a)
#' wktview(a, zoom=9)
#' }
wktview <- function(x, center = NULL, zoom = 5, fmt = 16) {
  UseMethod("wktview")
}

#' @export
wktview.character <- function(x, center = NULL, zoom = 5, fmt = 16) {
  not_some(x)
  make_view(x, center = center, zoom = zoom, fmt = fmt)
}

make_view <- function(x, center = NULL, zoom = 5, fmt = 16) {
  chek_for_pkg("leaflet")
  geojson <- wkt2geojson(x, fmt = fmt) %>% properties(style = list(NULL))
  cen <- centroid(geojson, center)
  leaflet::leaflet() %>%
    leaflet::addTiles() %>%
    leaflet::setView(lng = cen[1], lat = cen[2], zoom = zoom) %>%
    leaflet::addGeoJSON(geojson)
}

not_some <- function(x) {
  types <- c("POINT",'MULTIPOINT',"POLYGON","MULTIPOLYGON",
             "LINESTRING","MULTILINESTRING","GEOMETRYCOLLECTION")
  if (!grepl(paste0(types, collapse = "|"), x)) {
    stop("Only these supported:\n", paste0(types, collapse = ", "), call. = FALSE)
  }
}
