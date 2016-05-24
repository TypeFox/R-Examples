
#' Create a mesh of evenly spaced lines in another projection.
#'
#' @param x object to build line mesh for
#' @param proj the other projection
#'
#' @return spatial object
#' @export
#'
#' @examples
#' \dontrun{
#' library(maptools)
#' data(wrld_simpl)
#' library(raster)
#' w <- subset(wrld_simpl, NAME == "Australia")
#' plot(w)
#' laea <- pathologicule(w, "+proj=laea +lon_0=147 +lat_0=-42 +ellps=WGS84")
#' stere <- pathologicule(w, "+proj=stere +lon_0=147 +lat_0=-42 +ellps=WGS84")
#' plot(laea, add = TRUE, col = "dodgerblue")
#' plot(stere, add = TRUE, col = "firebrick")
#' pst90 <- "+proj=stere +lat_0=-90 +ellps=WGS84"
#' p <- spTransform(subset(wrld_simpl, coordinates(wrld_simpl)[,2] < -20), pst90)
#' plot(extent(p) + 1e6, asp = 1, type = "n"); plot(p, add = TRUE)
#' laea <- pathologicule(p, "+proj=laea +lon_0=147 +lat_0=-72 +ellps=WGS84")
#' stere <- pathologicule(p, "+proj=stere +lon_0=147 +lat_0=-42 +ellps=WGS84")
#' plot(laea, add = TRUE, col = "dodgerblue")
#' plot(stere, add = TRUE, col = "firebrick")
#' }
#' @importFrom raster projection projection<-
pathologicule <- function(x, proj) {
  y <- spTransform(x, proj)
  xr <- xrange(y)
  yr <- yrange(y)
  xs <- seq(xr[1], xr[2], length = 15)
  ys <- seq(yr[1], yr[2], length = 15)
  g <- graticule(xs, ys)
  projection(g) <- proj
  spTransform(g, projection(x))
}

