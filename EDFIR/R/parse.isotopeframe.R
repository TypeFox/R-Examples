parse.isotopeframe <-
function(isotope.frame) {
  ## parse the isotope frame produced by build.isotopeframe
  ## the isotope frame is a s x (2*d + 1) data frame, where s is the number of organisms
  ##    and d is the number of isotopes
  ## isotope.frame is organised like so: organism name, mean, sd, mean, sd, etc.
  ## NOTE: isotope input assumes a diagonal covariance matrix
  nrows = dim(isotope.frame)[1]
  ncols = dim(isotope.frame)[2]
  d = (ncols-1)/2
  mu = array(0, dim=c(nrows,d))
  sigma = array(0, dim=c(nrows, d))
  for (i in 1:nrows) {
    for (j in 1:d) {
      mu[i,j] = isotope.frame[i, 2*j]
      sigma[i,j] = isotope.frame[i, 1+2*j]
    }
  }
  list(mu=mu, sigma=sigma)
}
