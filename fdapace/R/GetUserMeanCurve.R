GetUserMeanCurve <- function (optns, obsGrid, regGrid, buff) {
  # If the user provided a mean function use it

  userMu = optns$userMu
  rangeUser <- range(optns$userMu$t)
  rangeObs <- range(obsGrid)
  if( rangeUser[1] > rangeObs[1] + buff || 
      rangeUser[2] < rangeObs[2] - buff   ) {
    stop('The range defined by the user provided mean does not cover the support of the data.')
  }

  mu = spline(userMu$t, userMu$mu, xout= obsGrid)$y
  muDense = spline(obsGrid,mu, xout=regGrid)$y
  bw_mu = NULL
 
  result <- list( 'mu' = mu, 'muDense'= muDense, 'bw_mu' = bw_mu)
  class(result) <- "SMC"
  return(result)
}
