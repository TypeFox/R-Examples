## special versions of approx:
##  approxTime:  interpolation of complete rows of a matrix or data frame
##  approxTime1: special case with one row only (slightly faster)

approxTime1 <- function (x, xout, rule = 1) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if ((!is.numeric(xout)) | (length(xout) != 1))
    stop("xout must be a scalar numeric value")
  if ((!is.numeric(rule)) | (length(rule) != 1))
    stop("rule must be a scalar numeric value")

  n <- nrow(x)
  if (xout >= x[n, 1]) {
    y <- c(xout, x[n, -1])
    if (rule == 1 & (xout > x[n + 1]))
      y[2:length(y)] <- NA
  }
  else if (xout <= x[1, 1]) {
    y <- c(xout, x[1, -1])
    if (rule == 1 & (xout < x[1]))
      y[2:length(y)] <- NA
  }
  else {
    i <- which.max(x[, 1] > xout)
    x1 <- x[i - 1, 1]
    x2 <- x[i, 1]
    y1 <- x[i - 1, ]
    y2 <- x[i, ]
    y <- y1 + (y2 - y1) * (xout - x1)/(x2 - x1)
  }
  names(y) <- dimnames(x)[[2]]
  y
}
