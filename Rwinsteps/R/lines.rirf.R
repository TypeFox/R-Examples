lines.rirf <- function(x, theta = x$theta, lwd = 2,
  col = "r", ...) {

  col <- getcolors(col, ncol(x$p))

  for(i in 1:ncol(x$p))
    lines.default(theta, x$p[, i], lwd = lwd, col = col[i], ...)
}
