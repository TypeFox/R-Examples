# Calculate a series of estimates of overlap based on bootstrap samples.

# Amat, Bmat: matrices with bootstrap samples in columns
# kmax: maximum value of k for optimal bandwidth estimation
# adjust: smoothing adjustment, a vector of length 3 for
#   Dhat1, Dhat4, Dhat5; if NA, the corresponding estimator is
#   ommitted. (adjust = 1/c in old code)
# n.grid: number of points at which to estimate density

# Returns: a matrix of overlap estimates with a column for
#   each estimator


bootEst <-
function(Amat, Bmat, kmax=3, adjust=c(0.8, 1, 4), n.grid=128) {
  nboot <- min(ncol(Amat), ncol(Bmat))
  bsamp <- matrix(NA, nboot, 3)
  colnames(bsamp) <- c("Dhat1", "Dhat4", "Dhat5")
  for(i in 1:nboot)
    bsamp[i, ] <- overlapEst(Amat[, i], Bmat[, i], kmax=kmax,
                          adjust=adjust, n.grid=n.grid) 
  return(bsamp)
}
