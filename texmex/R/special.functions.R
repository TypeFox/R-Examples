##' Miscellaneous internal functions
NULL

##' Accurately compute (exp(x) - 1) / x
##' @param x numeric vector
##' @return numeric vector
.exprel <- function(x) {
  .Call(.c.exprel, x)
}

##' Accurately compute log(1 + x) / x
##' @param x numeric vector
##' @return numeric vector
.log1prel <- function(x) {
  .Call(.c.log1prel, x)
}

##' Compute pmax(x y, -1) in such a way that zeros in x beat
##' infinities in y.
##'
##' This is a common pattern in much of the distribution code, so it's
##' worth factoring out.
##' @param x a numeric vector
##' @param y a numeric vector
##' @return an appropriate numeric vector
.specfun.safe.product <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  xy <- x * y
  xy[(x==0) & is.infinite(y)] <- 0
  pmax(xy, -1)
}


