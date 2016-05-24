#' To WKT
#'
#' @export
#' @param x Input
#' @examples
#' library("rgeos")
#' wkt <- "POLYGON((-180 -20, -140 55, 10 0, -140 -60, -180 -20))"
#' poly <- readWKT(wkt)
#' polys <- chop(x = poly)
#' to_wkt(polys)
#' to_wkt(polys)[[2]]
to_wkt <- function(x) {
  UseMethod("to_wkt")
}

#' @export
to_wkt.axewkt <- function(x) x

#' @export
to_wkt.SpatialPolygons <- function(x) {
  towkt(lapply(x@polygons, function(z) rgeos::writeWKT(sp::SpatialPolygons(list(z)))))
}

towkt <- function(x) {
  structure(x, class = "axewkt")
}
