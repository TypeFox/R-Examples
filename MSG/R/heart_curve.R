#' Draw a heart curve
#'
#' Calculate the coordinates of a heart shape and draw it with a polygon.
#' @param n the number of points to use when calculating the coordinates of the
#'   heart shape
#' @param ... other arguments to be passed to \code{\link[graphics]{polygon}},
#'   e.g. the color of the polygon (usually red)
#' @return NULL
#' @author Yihui Xie <\url{http://yihui.name}>
#' @export
#' @examples heart_curve()
#' heart_curve(col = 'red')
#' heart_curve(col = 'pink', border = 'red')
heart_curve = function(n = 101, ...) {
  y0 = seq(0, pi/2, length = n)
  x0 = sin(y0)
  y0 = c(y0, sqrt(1/4 - c(rev(x0) - .5)^2) + pi/2)
  x0 = c(x0, rev(x0))
  x0 = c(x0, -x0[c((2*n):(n+1), n:1)])
  y0 = c(y0, rev(y0))
  par(mar = rep(.05, 4))
  plot(x0, y0, type = 'n', ann = FALSE, axes = FALSE)
  polygon(x0, y0, ...)
}
