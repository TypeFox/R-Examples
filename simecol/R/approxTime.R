## special versions of approx:
##  approxTime:  interpolation of complete rows of a matrix or data frame
##  approxTime1: special case with one row only (slightly faster)

approxTime <- function(x, xout, ...) {
  if (is.data.frame(x)) {x <- as.matrix(x); wasdf <- TRUE} else wasdf <- FALSE
  if (!is.matrix(x)) stop("x must be a matrix or data frame")
  m <- ncol(x)
  y <- matrix(0, nrow=length(xout), ncol=m)
  y[,1] <- xout
  for (i in 2:m) {
    y[,i] <- as.vector(approx(x[,1], x[,i], xout, ...)$y)
  }
  if (wasdf) y <- as.data.frame(y)
  names(y) <- dimnames(x)[[2]]
  y
}
