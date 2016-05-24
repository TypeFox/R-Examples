whiteKernDiagCompute <-
function (kern, x) {
  k <- matrix(kern$variance, dim(as.array(x))[1], 1)
  return (k)
}
