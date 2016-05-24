lines.riif <- function(x, theta = x$theta, lwd = 2,
  col = "r", ...) {

  col <- getcolors(col, ncol(x$pq))

  for(i in 1:ncol(x$pq))
    lines.default(theta, x$pq[, i], lwd = lwd, col = col[i], ...)
}
