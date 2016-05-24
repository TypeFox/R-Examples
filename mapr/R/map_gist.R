#' Make an interactive map to view in the browser as a GitHub gist
#'
#' @export
#' @template args
#' @param description Description for the Github gist, or leave to default (=no description)
#' @param public (logical) Whether gist is public (default: TRUE)
#' @param browse If TRUE (default) the map opens in your default browser.
#' @param ... Further arguments passed on to \code{\link{style_geojson}}
#'
#' @details See \code{\link[gistr]{gist_auth}} for help on authentication
#'
#' @examples \dontrun{
#' ## spocc
#' library("spocc")
#' spp <- c('Danaus plexippus', 'Accipiter striatus', 'Pinus contorta')
#' dat <- occ(spp, from=c('gbif','ecoengine'), limit=30, gbifopts=list(hasCoordinate=TRUE))
#' dat <- fixnames(dat, "query")
#'
#' # Define colors
#' map_gist(dat, color=c('#976AAE','#6B944D','#BD5945'))
#' map_gist(dat$gbif, color=c('#976AAE','#6B944D','#BD5945'))
#' map_gist(dat$ecoengine, color=c('#976AAE','#6B944D','#BD5945'))
#'
#' # Define colors and marker size
#' map_gist(dat, color=c('#976AAE','#6B944D','#BD5945'), size=c('small','medium','large'))
#'
#' # Define symbols
#' map_gist(dat, symbol=c('park','zoo','garden'))
#'
#' ## rgbif
#' library("rgbif")
#' res <- occ_search(scientificName = "Puma concolor", limit = 100)
#' map_gist(res)
#'
#' ## data.frame
#' df <- data.frame(name = c('Poa annua', 'Puma concolor', 'Foo bar'),
#'                  longitude = c(-120, -121, -121),
#'                  latitude = c(41, 42, 45), stringsAsFactors = FALSE)
#' map_gist(df)
#'
#' ### usage of occ2sp()
#' #### SpatialPoints
#' spdat <- occ2sp(dat)
#' map_gist(spdat)
#' #### SpatialPointsDataFrame
#' spdatdf <- as(spdat, "SpatialPointsDataFrame")
#' map_gist(spdatdf)
#' }
map_gist <- function(x, description = "", public = TRUE, browse = TRUE,
                     lon = 'longitude', lat = 'latitude', ...) {
  UseMethod("map_gist")
}

#' @export
map_gist.occdat <- function(x, description = "", public = TRUE, browse = TRUE,
                            lon = 'longitude', lat = 'latitude', ...) {
  map_gister(spocc::occ2df(x), description, public, browse, ...)
}

#' @export
map_gist.occdatind <- function(x, description = "", public = TRUE, browse = TRUE,
                               lon = 'longitude', lat = 'latitude', ...) {
  map_gister(spocc::occ2df(x), description, public, browse, ...)
}

#' @export
map_gist.gbif <- function(x, description = "", public = TRUE, browse = TRUE,
                          lon = 'longitude', lat = 'latitude', ...) {
  x <- x$data
  x <- re_name(x, c('decimalLatitude' = 'latitude'))
  x <- re_name(x, c('decimalLongitude' = 'longitude'))
  map_gister(x, description, public, browse, ...)
}

#' @export
map_gist.data.frame <- function(x, description = "", public = TRUE, browse = TRUE,
                                lon = 'longitude', lat = 'latitude', ...) {

  x <- guess_latlon(x, lat, lon)
  map_gister(x, description, public, browse, ...)
}

#' @export
map_gist.SpatialPoints <- function(x, description = "", public = TRUE, browse = TRUE,
                                   lon = 'longitude', lat = 'latitude', ...) {
  x <- data.frame(x)
  x <- guess_latlon(x, lat, lon)
  map_gister(x, description, public, browse, ...)
}

#' @export
map_gist.SpatialPointsDataFrame <- function(x, description = "", public = TRUE, browse = TRUE,
                                            lon = 'longitude', lat = 'latitude', ...) {
  x <- data.frame(x)
  x <- guess_latlon(x, lat, lon)
  map_gister(x, description, public, browse, ...)
}

#' @export
map_gist.default <- function(x, description = "", public = TRUE, browse = TRUE,
                             lon = 'longitude', lat = 'latitude', ...) {
  stop(sprintf("map_gist does not support input of class '%s'", class(x)), call. = FALSE)
}

# helpers
map_gister <- function(x, description, public, browse, ...) {
  datgeojson <- style_geojson(input = x, var = "name", ...)
  file <- tempfile(fileext = ".csv")
  write.csv(datgeojson, file)
  geofile <- togeojson2(file)
  gistr::gist_create(geofile, description = description, public = public, browse = browse)
}
