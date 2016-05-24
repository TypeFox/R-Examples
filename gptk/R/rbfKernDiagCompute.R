rbfKernDiagCompute <-
function (kern, x) {
  k <- matrix(kern$variance, dim(as.array(x))[1], 1)

  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    k <- k * sqrt(kern$inverseWidth/(2*pi))

  return (k)
}
