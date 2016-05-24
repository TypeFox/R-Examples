
#' Overlay points on an OSM plot
#'
#' Plot points on a plot created by \link{osm.plot}. This is a simple wrapper around
#' \code{points()}.
#'
#' @param x X coordinate vector or object as parsed by \code{xy.coords}
#' @param y Y coordinate vector
#' @param epsg EPSG code of the supplied coordinates
#' @param toepsg EPSG code of the projected coordinates to be ploted
#' @param ... Args passed on to \code{points}
#'
#' @export
#'
#' @examples
#' library(rosm)
#' library(prettymapr)
#' locs <- geocode(c("wolfville, ns", "kentville, ns", "halifax, ns"))
#' prettymap({
#'   osm.plot(searchbbox("nova scotia"))
#'   osm.points(locs$lon, locs$lat, pch=18, cex=0.7)
#' })
#'
#'
osm.points <- function(x, y=NULL, epsg=4326, toepsg=3857, ...) {
  graphics::points(.projpts(x, y, epsg, toepsg), ...)
}

#' Overlay segments on an OSM plot
#'
#' Plot segments on a plot created by \link{osm.plot}. This is a simple wrapper around
#' \code{segments()}.
#'
#' @param x0 X1 coordinate vector
#' @param y0 Y1 coordinate vector
#' @param x1 X2 coordinate vector
#' @param y1 Y2 coordinate vector
#' @param epsg EPSG code of the supplied coordinates
#' @param toepsg EPSG code of the projected coordinates to be ploted
#' @param ... Args passed on to \code{points}
#'
#' @export
#'
#' @examples
#' library(rosm)
#' library(prettymapr)
#' locs <- geocode(c("wolfville, ns", "kentville, ns", "halifax, ns"))
#' prettymap({
#'   osm.plot(searchbbox("nova scotia"))
#'   osm.segments(locs$lon[1:2], locs$lat[1:2], locs$lon[2:3], locs$lat[2:3])
#' })
#'
#'
osm.segments <- function(x0, y0, x1=x0, y1=y0, epsg=4326, toepsg=3857, ...) {
  c0 <- .projpts(x0, y0, epsg, toepsg)
  c1 <- .projpts(x1, y1, epsg, toepsg)
  graphics::segments(c0[,1], c0[,2], c1[,1], c1[,2])
}

#' Overlay lines on an OSM plot
#'
#' Plot lines on a plot created by \link{osm.plot}. This is a simple wrapper around
#' \code{points()}.
#'
#' @param x X coordinate vector or object as parsed by \code{xy.coords}
#' @param y Y coordinate vector
#' @param epsg EPSG code of the supplied coordinates
#' @param toepsg EPSG code of the projected coordinates to be ploted
#' @param ... Args passed on to \code{lines}
#'
#' @export
#'
#' @examples
#' library(rosm)
#' library(prettymapr)
#' locs <- geocode(c("wolfville, ns", "kentville, ns", "halifax, ns"))
#' prettymap({
#'   osm.plot(searchbbox("nova scotia"))
#'   osm.lines(locs$lon, locs$lat, lwd=2)
#' })
#'
osm.lines <- function(x, y=NULL, epsg=4326, toepsg=3857, ...) {
  graphics::lines(.projpts(x, y, epsg, toepsg), ...)
}

#' Overlay a polygon on an OSM plot
#'
#' Plot a polygon on a plot created by \link{osm.plot}. This is a simple wrapper around
#' \code{polygon()}.
#'
#' @param x X coordinate vector or object as parsed by \code{xy.coords}
#' @param y Y coordinate vector
#' @param epsg EPSG code of the supplied coordinates
#' @param toepsg EPSG code of the projected coordinates to be ploted
#' @param ... Args passed on to \code{polygon}
#'
#' @export
#'
#' @examples
#' library(rosm)
#' library(prettymapr)
#' locs <- geocode(c("wolfville, ns", "kentville, ns", "halifax, ns"))
#' prettymap({
#'   osm.plot(searchbbox("nova scotia"))
#'   osm.polygon(locs$lon, locs$lat)
#' })
#'
osm.polygon <- function(x, y=NULL, epsg=4326, toepsg=3857, ...) {
  graphics::polygon(.projpts(x, y, epsg, toepsg), ...)
}


.projpts <- function(x, y=NULL, epsg, toepsg) {
  rgdal::CRSargs(sp::CRS("+init=epsg:3857")) #hack to load rgdal namespace
  xy <- grDevices::xy.coords(x, y)
  sp::coordinates(sp::spTransform(sp::SpatialPoints(sp::coordinates(cbind(xy$x, xy$y)),
                                    proj4string=sp::CRS(sprintf("+init=epsg:%s", epsg))),
                  sp::CRS(sprintf("+init=epsg:%s", toepsg))))
}

