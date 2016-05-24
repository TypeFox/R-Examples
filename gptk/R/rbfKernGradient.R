rbfKernGradient <-
function (kern, x, x2, covGrad) {
  if ( nargs()==3 ) {
    k <- rbfKernCompute(kern, x)
    dist2xx <- .dist2(x, x)
    covGrad <- x2
  } else if ( nargs()==4 ) {
    k <- rbfKernCompute(kern, x, x2)
    dist2xx <- .dist2(x, x2)
  }

  g <- array()
  if ("isNormalised" %in% names(kern) && kern$isNormalised) {
    g[1] <- -0.5*sum(covGrad*k*dist2xx) +
      0.5 * sum(covGrad*k)/kern$inverseWidth
  }
  else {
    g[1] <- -0.5*sum(covGrad*k*dist2xx)
  }
  g[2] <- sum(covGrad*k)/kern$variance

  if ( any(is.nan(g)) )
    warning("g is NaN.")

  return (g)
}
