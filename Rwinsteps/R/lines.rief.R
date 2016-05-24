lines.rief <- function(x, theta = x$theta, lwd = 2,
  col = "r", ...) {

  col <- getcolors(col, ncol(x$e))

  for(i in 1:ncol(x$e))
    lines.default(theta, x$e[, i], lwd = lwd, col = col[i], ...)
}
