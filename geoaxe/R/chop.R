#' Split polygon into many
#'
#' @export
#' @param x Spatial object
#' @param size size of each side of each cell, which makes a square cell
#' @param n number of cells to make in each dimension, same number used for
#' each dimension
#' @details Works on spatial classes of type \code{SpatialPolygons},
#' Well-Known Text character strings, and GeoJSON character strings and lists
#' @examples
#' library("rgeos")
#' wkt <- "POLYGON((-180 -20, -140 55, 10 0, -140 -60, -180 -20))"
#'
#' # SpatialPolygons input
#' poly <- readWKT(wkt)
#' plot(poly)
#' polys <- chop(x = poly)
#' to_wkt(polys)
#' to_wkt(polys)[[2]]
#' plot(polys)
#' plot(poly, add = TRUE, lwd = 6)
#'
#' # SpatialPolygonsDataFrame input
#' class(poly)
#' polydf <- as(poly, "SpatialPolygonsDataFrame")
#' class(polydf)
#' chop(polydf)
#'
#' # WKT character input
#' chop(wkt)
#'
#' # geojson character input
#' file <- system.file("examples", "sample1.geojson", package = "geoaxe")
#' x <- readLines(file)
#' chop(x)
#'
#' # geojson json input
#' x <- structure(x, class = "json")
#' chop(x)
chop <- function(x, size = 10, n = 20) {
  UseMethod("chop")
}

#' @export
chop.SpatialPolygons <- function(x, size = 10, n = 20) {
  box <- sp::bbox(x)
  gt <- sp::GridTopology(c(box[1,1], box[2,1]), rep(size, 2), rep(n, 2))
  gr <- as(as(sp::SpatialGrid(gt), "SpatialPixels"), "SpatialPolygons")
  rgeos::gIntersection(x, gr, byid = TRUE, drop_lower_td = TRUE)
}

#' @export
chop.SpatialPolygonsDataFrame <- function(x, size = 10, n = 20) {
  chop(as(x, "SpatialPolygons"))
}

#' @export
chop.character <- function(x, size = 10, n = 20) {
  switch(wkt_geojson(x),
    wkt = chop(rgeos::readWKT(x), size = size, n = n),
    geojson = as_SpatialPolygons(jsonlite::fromJSON(x, FALSE))
  )
}

#' @export
chop.json <- function(x, size = 10, n = 20) {
  chop(unclass(x))
}

#' @export
chop.list <- function(x, size = 10, n = 20) {
  ## FIXME, need to check that list is geojson type
  as_SpatialPolygons(x)
}

#' @export
chop.default <- function(x, size = 10, n = 20) {
  stop(sprintf("`chop()` method not implemented for %s.", class(x)), call. = FALSE)
}

wkt_geojson <- function(x) {
  if (grepl("\\{|\\}", x)) {
    "geojson"
  } else if (grepl("POLYGON|MULTIPOLYGON", x)) {
    "wkt"
  } else {
    stop("input must be WKT or GeoJSON", call. = FALSE)
  }
}
