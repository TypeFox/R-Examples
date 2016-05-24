rbfKernCompute <-
function (kern, x, x2=NULL) {
  if ( nargs() < 3 ) {
    n2 <- .dist2(x,x)
  } else {
    n2 <- .dist2(x,x2)
  }

  wi2 <- 0.5*kern$inverseWidth
  k <- kern$variance*exp(-n2*wi2)

  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    k <- k * sqrt(kern$inverseWidth/(2*pi))
  
  return (k)
}
