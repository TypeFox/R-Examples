

rvngamma <- function (n=1, shape, rate = 1, scale = 1/rate) {
  if (any(shape < 0, na.rm=TRUE)) {
    stop("shape parameter must be nonnegative")
  }
  shape <- (1/3 + shape)
  rvvapply(stats:::rgamma, n.=n, shape=shape, scale=scale)
}

