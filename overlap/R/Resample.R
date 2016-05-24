
# Nonparametric resampling

# x: vector of observation times in radians
# nb: number of bootstrap samples required
# smooth: if TRUE, smoothed resamples are returned
# kmax, adjust, n.grid: arguments passed to densityFit

resample <-
function(x, nb, smooth=TRUE, kmax=3, adjust=1, n.grid=512)  {
  n <- length(x)
  if(smooth) {
    bw0 <- getBandWidth(x, kmax=kmax)
    if(is.na(bw0))  # if Uniroot Error
      return(NA)
    dA <- densityFit(x, seq(0, 2*pi, length=n.grid), bw0 / adjust)
    res <- matrix(rejectSampleRad(n * nb, dA), n, nb)
  } else {
    res <- matrix(sample(x, n * nb, replace=TRUE), n, nb)
  }
  return(res)
}


rejectSampleRad <-
function(n, dA) {
  # Uses rejection sampling to draw a sample from a circular density. 
  # Args:
  #   n = sample size required
  #   dA = a density on a regular grid from 0 to 2*pi
  # Returns:
  #   A sample length n in [0, 2*pi)
  #
  uplim <- max(dA)
  nx <- length(dA) - 1    # First & last entries both refer to midnight
  iCand <- sample.int(nx, 2 * n, replace=TRUE)  # these are INDICES
  accept <- dA[iCand] > runif(2 * n, 0, uplim)
  iSamp <- iCand[accept]
  while(length(iSamp) < n)  { 
    iCand <- sample.int(nx, n, replace=TRUE) 
    accept <- dA[iCand] > runif(n, 0, uplim)
    iSamp <- c(iSamp, iCand[accept])
  }
  grid <- seq(0, 2*pi, length=length(dA))
  return(grid[iSamp[1:n]])  
}
