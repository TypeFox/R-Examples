#' @importFrom raster xmin xmax
xrange <- function(x) {
  c(xmin(x), xmax(x))
}
#' @importFrom raster ymin ymax
yrange <- function(x) {
  c(ymin(x), ymax(x))
}
