#' ggplot2 mapping
#'
#' @export
#' @template args
#' @param map (character) One of world, world2, state, usa, county, france, italy, or nz
#' @param point_color Default color of your points
#' @param ... Ignored
#' @return A ggplot2 map, of class \code{gg/ggplot}
#' @examples \dontrun{
#' ## spocc
#' library("spocc")
#' dat <- occ(query = 'Lynx rufus californicus', from = 'gbif', limit=100)
#' map_ggplot(dat)
#' map_ggplot(dat$gbif)
#' map_ggplot(dat, "usa")
#' map_ggplot(dat, "county")
#'
#' ### usage of occ2sp()
#' #### SpatialPoints
#' spdat <- occ2sp(dat)
#' map_ggplot(spdat)
#' #### SpatialPointsDataFrame
#' spdatdf <- as(spdat, "SpatialPointsDataFrame")
#' map_ggplot(spdatdf)
#'
#' ## rgbif
#' library("rgbif")
#' res <- occ_search(scientificName = "Puma concolor", limit = 100)
#' map_ggplot(res)
#'
#' ## data.frame
#' df <- data.frame(name = c('Poa annua', 'Puma concolor', 'Foo bar'),
#'                  longitude = c(-120, -121, -121),
#'                  latitude = c(41, 42, 45), stringsAsFactors = FALSE)
#' map_ggplot(df)
#'}
map_ggplot <- function(x, map = "world", point_color = "#86161f",
                       lon = 'longitude', lat = 'latitude', ...) {
  UseMethod("map_ggplot")
}

#' @export
map_ggplot.occdat <- function(x, map = "world", point_color = "#86161f",
                              lon = 'longitude', lat = 'latitude', ...) {

  x <- spocc::occ2df(x)
  make_amap(dat_cleaner(x, lon = 'longitude', lat = 'latitude'), map, point_color)
}

#' @export
map_ggplot.occdatind <- function(x, map = "world", point_color = "#86161f",
                              lon = 'longitude', lat = 'latitude', ...) {
  x <- spocc::occ2df(x)
  make_amap(dat_cleaner(x, lon = 'longitude', lat = 'latitude'), map, point_color)
}

#' @export
map_ggplot.gbif <- function(x, map = "world", point_color = "#86161f",
                            lon = 'longitude', lat = 'latitude', ...) {
  make_amap(dat_cleaner(x$data, lon = 'decimalLongitude', lat = 'decimalLatitude'), map, point_color)
}

#' @export
map_ggplot.SpatialPoints <- function(x, map = "world", point_color = "#86161f",
                    lon = 'longitude', lat = 'latitude', ...) {
  make_amap(data.frame(x), map, point_color)
}

#' @export
map_ggplot.SpatialPointsDataFrame <- function(x, map = "world", point_color = "#86161f",
                    lon = 'longitude', lat = 'latitude', ...) {
  make_amap(data.frame(x), map, point_color)
}

#' @export
map_ggplot.data.frame <- function(x, map = "world", point_color = "#86161f",
                                  lon = 'longitude', lat = 'latitude', ...) {
  make_amap(dat_cleaner(x, lon = 'longitude', lat = 'latitude'), map, point_color)
}

#' @export
map_ggplot.default <- function(x, map = "world", point_color = "#86161f",
                               lon = 'longitude', lat = 'latitude', ...) {
  stop(sprintf("map_ggplot does not support input of class '%s'", class(x)), call. = FALSE)
}

### helpers ------------------------------------------
sutils_blank_theme <- function(){
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        plot.margin = rep(ggplot2::unit(0, "null"), 4))
}

make_amap <- function(x, map, point_color) {
  wmap <- suppressMessages(ggplot2::map_data(map))
  latitude <- longitude <- lat <- long <- decimalLongitude <- decimalLatitude <- group <- NA
  ggplot(x, aes(longitude, latitude)) +
    geom_point(color = point_color, size = 3) +
    geom_polygon(aes(long, lat, group = group), fill = NA, colour = "black", data = wmap) +
    sutils_blank_theme()
}
