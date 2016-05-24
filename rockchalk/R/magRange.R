##' magRange
##' Magnify the range of a variable.
##'
##' By default, R's range function returns the minimum and maximum
##' values of a variable. This returns a magnified range. It is used
##' for some plotting functions in the rockchalk package
##'
##' @usage magRange(x, mult = 1.25)
##' @param x an R vector variable
##' @param mult a multiplier by which to magnify the range of the
##' variable.  A value of 1 leaves the range unchanged. May be a
##' scalar, in which case both ends of the range are magnified by the
##' same amount.  May also be a two valued vector, such as c(minMag,
##' maxMag), in which case the magnification applied to the minimum is
##' minMag and the magnification of the maximum is maxMag.  After
##' version 1.5.5, mult may be smaller than 1, thus shrinking the
##' range. Setting mult to values closer to 0 causes the range to
##' shrink to the center point from both sides.
##' @name magRange
##' @export magRange
##' @examples
##' x1 <- c(0, 0.5, 1.0)
##' range(x1)
##' magRange(x1, 1.1)
##' magRange(x1, c(1.1, 1.4))
##' magRange(x1, 0.5)
##' magRange(x1, c(0.1, 0.1))
##' x1 <- rnorm(100)
##' range(x1)
##' magRange(x1)
##' magRange(x1, 1.5)
##' magRange(x1, c(1,1.5))

magRange <- function(x, mult = 1.25){
  xr <- range(x, na.rm = TRUE)
  xm <- 0.5*(xr[1] + xr[2])
  magXr <- vector(mode = "numeric", length = 2)

  if (length(mult) > 2) stop("magRange: supply only 2 multipliers")
  if (length(mult) == 1) mult <- c(mult, mult)

  if ( mult[2] > 1 ) {
      magXr[2] <- xr[2] + 2 * (mult[2] - 1)*(xr[2] - xm)
  } else {
      magXr[2] <- xm + 0.5 * mult[2] * (xr[2] -xr[1])
  }

  if ( mult[1] > 1 ) {
      magXr[1] <- xr[1] - 2 * (mult[1] - 1) *( xm - xr[1])
  } else {
      magXr[1] <- xm - 0.5 * mult[1] * (xr[2] -xr[1])
  }
  magXr
}


