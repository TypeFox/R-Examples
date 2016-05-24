#' ggpmap visualization of species occurences
#'
#' @export
#' @template args
#' @param zoom zoom level for map. Adjust depending on how your data look.
#' @param point_color Default color of your points
#' @param ... Ignored
#' @examples \dontrun{
#' ## spocc
#' library("spocc")
#' gd <- occ(query = 'Accipiter striatus', from = 'gbif', limit=75, has_coords = TRUE)
#' map_ggmap(gd)
#' map_ggmap(gd$gbif)
#'
#' ## rgbif
#' library("rgbif")
#' res <- occ_search(scientificName = "Puma concolor", limit = 100)
#' map_ggmap(res)
#'
#' ## data.frame
#' df <- data.frame(name = c('Poa annua', 'Puma concolor', 'Foo bar'),
#'                  longitude = c(-120, -121, -123),
#'                  latitude = c(41, 42, 45), stringsAsFactors = FALSE)
#' map_ggmap(df)
#'
#' ### usage of occ2sp()
#' #### SpatialPoints
#' spdat <- occ2sp(gd)
#' map_ggmap(spdat)
#' #### SpatialPointsDataFrame
#' spdatdf <- as(spdat, "SpatialPointsDataFrame")
#' map_ggmap(spdatdf)
#'}
map_ggmap <- function(x, zoom = 3, point_color = "#86161f",
                      lon = 'longitude', lat = 'latitude', ...) {
  UseMethod("map_ggmap")
}

#' @export
map_ggmap.occdat <- function(x, zoom = 3, point_color = "#86161f",
                             lon = 'longitude', lat = 'latitude', ...) {
  x <- spocc::occ2df(x)
  map_ggmapper(x, zoom, point_color)
}

#' @export
map_ggmap.occdatind <- function(x, zoom = 3, point_color = "#86161f",
                                lon = 'longitude', lat = 'latitude', ...) {
  x <- spocc::occ2df(x)
  map_ggmapper(x, zoom, point_color)
}

#' @export
map_ggmap.gbif <- function(x, zoom = 3, point_color = "#86161f",
                           lon = 'longitude', lat = 'latitude', ...) {
  x <- guess_latlon(x$data, lon = 'decimalLongitude', lat = 'decimalLatitude')
  map_ggmapper(x, zoom, point_color)
}

#' @export
map_ggmap.data.frame <- function(x, zoom = 3, point_color = "#86161f",
                                 lon = 'longitude', lat = 'latitude', ...) {
  x <- guess_latlon(x, lat, lon)
  map_ggmapper(x, zoom, point_color)
}

#' @export
map_ggmap.SpatialPoints <- function(x, zoom = 3, point_color = "#86161f",
                                    lon = 'longitude', lat = 'latitude', ...) {
  x <- data.frame(x)
  x <- guess_latlon(x, lat, lon)
  map_ggmapper(x, zoom, point_color)
}

#' @export
map_ggmap.SpatialPointsDataFrame <- function(x, zoom = 3, point_color = "#86161f",
                                             lon = 'longitude', lat = 'latitude', ...) {
  x <- data.frame(x)
  x <- guess_latlon(x, lat, lon)
  map_ggmapper(x, zoom, point_color)
}

#' @export
map_ggmap.default <- function(x, zoom = 3, point_color = "#86161f",
                              lon = 'longitude', lat = 'latitude', ...) {
  stop(sprintf("map_ggmap does not support input of class '%s'", class(x)), call. = FALSE)
}

## helpers ---------------------
map_center <- function(x) {
  min_lat <- min(x$latitude, na.rm = TRUE)
  max_lat <- max(x$latitude, na.rm = TRUE)
  min_long <- min(x$longitude, na.rm = TRUE)
  max_long <- max(x$longitude, na.rm = TRUE)
  center_lat <- min_lat + (max_lat - min_lat)/2
  center_long <- min_long + (max_long - min_long)/2
  c(lon = center_long, lat = center_lat)
}

map_ggmapper <- function(x, zoom, point_color) {
  check4pkg("ggmap")
  x <- x[complete.cases(x$latitude, x$longitude), ]
  x <- x[!x$latitude == 0 & !x$longitude == 0, ]
  species_map <- ggmap::get_map(location = map_center(x), zoom = zoom, maptype = "terrain")
  latitude <- longitude <- NA
  ggmap::ggmap(species_map) +
    ggplot2::geom_point(data = x[, c("latitude", "longitude")],
                        ggplot2::aes(x = longitude, y = latitude), color = point_color, size = 3) +
    ggplot2::ggtitle(paste0("Distribution of ", unique(x$name))) +
    ggplot2::labs(x = "Longitude", y = "Latitude")
}
