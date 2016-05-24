#' Make an interactive map with Leaflet.js
#'
#' @export
#'
#' @template args
#' @param ... Ignored
#' @examples \dontrun{
#' ## spocc
#' library("spocc")
#' (out <- occ(query='Accipiter striatus', from='gbif', limit=50, has_coords=TRUE))
#' ### with class occdat
#' map_leaflet(out)
#' ### with class occdatind
#' map_leaflet(out$gbif)
#' ### use occ2sp
#' map_leaflet(occ2sp(out))
#'
#' ## rgbif
#' library("rgbif")
#' res <- occ_search(scientificName = "Puma concolor", limit = 100)
#' map_leaflet(res)
#'
#' ## SpatialPoints class
#' df <- data.frame(longitude = c(-120,-121),
#'                  latitude = c(41, 42), stringsAsFactors = FALSE)
#' x <- SpatialPoints(df)
#' map_leaflet(x)
#'
#' ## SpatialPointsDataFrame class
#' library("rgbif")
#' res <- occ_search(scientificName = "Puma concolor", limit = 100)
#' x <- res$data
#' library("sp")
#' x <- x[complete.cases(x$decimalLatitude, x$decimalLongitude), ]
#' coordinates(x) <- ~decimalLongitude+decimalLatitude
#' map_leaflet(x)
#'
#' ## data.frame
#' df <- data.frame(name = c('Poa annua', 'Puma concolor'),
#'                  longitude = c(-120,-121),
#'                  latitude = c(41, 42), stringsAsFactors = FALSE)
#' map_leaflet(df)
#' }
map_leaflet <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  UseMethod("map_leaflet")
}

#' @export
map_leaflet.occdat <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  make_map_ll(dat_cleaner(spocc::occ2df(x), lon, lat))
}

#' @export
map_leaflet.occdatind <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  make_map_ll(dat_cleaner(spocc::occ2df(x), lon, lat))
}

#' @export
map_leaflet.SpatialPoints <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  make_map(x)
}

#' @export
map_leaflet.SpatialPointsDataFrame <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  make_map(x)
}

#' @export
map_leaflet.gbif <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  make_map_ll(dat_cleaner(x$data, lon = 'decimalLongitude', lat = 'decimalLatitude'))
}

#' @export
map_leaflet.data.frame <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  make_map_ll(dat_cleaner(x, lon, lat))
}

#' @export
map_leaflet.default <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  stop(sprintf("map_leaflet does not support input of class '%s'", class(x)), call. = FALSE)
}

# helpers ------------------------------------
dat_cleaner <- function(x, lon = 'longitude', lat = 'latitude') {
  x <- guess_latlon(x, lat, lon)
  x[complete.cases(x$latitude, x$longitude), ]
}

make_map <- function(x) {
  lf <- leaflet::leaflet(data = x)
  lf <- leaflet::addTiles(lf)
  leaflet::addMarkers(lf)
}

make_map_ll <- function(x) {
  lf <- leaflet::leaflet(data = x)
  lf <- leaflet::addTiles(lf)
  leaflet::addMarkers(lf, ~longitude, ~latitude)
}
