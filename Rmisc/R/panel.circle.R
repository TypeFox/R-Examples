#' Circle Drawing
#' 
#' A panel function for drawing circles.
#' 
#' @param x The x coordinate of the circle center
#' @param y The y coordinate of the circle center
#' @param r The radius of the circle
#' @param segments The number of polygon segments used to create the circle
#' @param groups A factor defining groups
#' @param ... Additional arguments passed to panel.polygon
#' 
#' @importFrom lattice panel.polygon panel.superpose
#' 
#' @export
#' 
#' @examples
#' panel.circle(0, 0, 10)
#' 
panel.circle <- function(x, y, r, segments = 50L, groups = NULL, ...)
{
  x <- as.numeric(x)
  y <- as.numeric(y)
  r <- as.numeric(r)
  if (!is.null(groups)) {
    panel.superpose(x = x, y = y, r = r, groups = groups, segments = segments,
                    panel.groups = panel.circle, ...)
  } else {
    angle.inc <- 2 * pi/segments
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    xv <- cos(angles) * r + x
    yv <- sin(angles) * r + y
    panel.polygon(xv, yv, ...)
  }
}