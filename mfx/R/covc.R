covc <-
function(cluster, estfun, N1, K1, NROW1, fm1) {
  cluster = factor(cluster)
  # Calculate the "meat" of the sandwich estimators
  # call this middle to avoid partial matching
  u = apply(estfun, 2, function(x) tapply(x, cluster, sum))
  middle = crossprod(u)/N1
  
  # Calculations for degrees-of-freedom corrections, followed 
  # by calculation of the variance-covariance estimate.
  # NOTE: NROW/N is a kluge to address the fact that sandwich uses the
  # wrong number of rows (includes rows omitted from the regression).
  M = length(levels(cluster))
  dfc = (M/(M-1)) * ((N1-1)/(N1-K1))
  return(dfc * NROW1/N1 * sandwich(fm1, meat.=middle))
}
