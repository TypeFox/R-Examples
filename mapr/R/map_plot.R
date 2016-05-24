#' Base R mapping
#'
#' @export
#' @template args
#' @param ... Further args to \code{\link{points}}
#' @examples \dontrun{
#' ## spocc
#' library("spocc")
#' (out <- occ(query='Accipiter striatus', from='gbif', limit=25, has_coords=TRUE))
#' ### class occdat
#' map_plot(out)
#' ### class occdatind
#' map_plot(out$gbif)
#'
#' ## rgbif
#' library("rgbif")
#' res <- occ_search(scientificName = "Puma concolor", limit = 100)
#' map_plot(res)
#'
#' ## data.frame
#' df <- data.frame(name = c('Poa annua', 'Puma concolor', 'Foo bar'),
#'                  longitude = c(-120, -121, -121),
#'                  latitude = c(41, 42, 45), stringsAsFactors = FALSE)
#' map_plot(df)
#'
#' ### usage of occ2sp()
#' #### SpatialPoints
#' spdat <- occ2sp(out)
#' map_plot(spdat)
#' #### SpatialPointsDataFrame
#' spdatdf <- as(spdat, "SpatialPointsDataFrame")
#' map_plot(spdatdf)
#' }
map_plot <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  UseMethod("map_plot")
}

#' @export
map_plot.occdat <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  df <- spocc::occ2df(x)
  plot_er(plot_prep(df), ...)
}

#' @export
map_plot.occdatind <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  df <- spocc::occ2df(x)
  plot_er(plot_prep(df), ...)
}

#' @export
map_plot.gbif <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  df <- x$data
  df <- df[complete.cases(df$decimalLatitude, df$decimalLongitude), ]
  df <- df[df$decimalLongitude != 0, ]
  sp::coordinates(df) <- ~decimalLongitude + decimalLatitude
  plot_er(df, ...)
}

#' @export
map_plot.data.frame <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  x <- guess_latlon(x, lat, lon)
  plot_er(plot_prep(x), ...)
}

#' @export
map_plot.SpatialPoints <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  x <- data.frame(x)
  x <- guess_latlon(x, lat, lon)
  plot_er(plot_prep(x), ...)
}

#' @export
map_plot.SpatialPointsDataFrame <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  x <- data.frame(x)
  x <- guess_latlon(x, lat, lon)
  plot_er(plot_prep(x), ...)
}

#' @export
map_plot.default <- function(x, lon = 'longitude', lat = 'latitude', ...) {
  stop(sprintf("map_plot does not support input of class '%s'", class(x)), call. = FALSE)
}

plot_prep <- function(x) {
  x <- x[complete.cases(x$latitude, x$longitude), ]
  x <- x[x$longitude != 0, ]
  sp::coordinates(x) <- ~longitude + latitude
  x
}

plot_er <- function(x, ...) {
  sp::proj4string(x) <- sp::CRS("+init=epsg:4326")
  sp::plot(rworldmap::getMap())
  graphics::points(x, col = "red", ...)
}
