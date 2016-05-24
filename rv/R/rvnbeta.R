

rvnbeta <- function (n=1, shape1, shape2) {
  if (any(shape1<0, na.rm=TRUE) || any(shape2<0, na.rm=TRUE)) {
    stop("Neutral Beta distribution requires positive parameters")
  }
  shape1 <- (shape1 + 1/3)
  shape2 <- (shape2 + 1/3)
  rvvapply(stats:::rbeta, n.=n, shape1=shape1, shape2=shape2)
}

