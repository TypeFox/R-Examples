getsphericalharmonicseven <- function( order, theta, phi) {
  
  ## compute spherical harmonics 
  ## (symmetric modified SH Basis used by Descoteaux (2008))
  
  order <- as.integer( max( 0, order))
  if ( order%%2 == 1) {
    warning("getsphericalharmonicseven: maximum order needs to be even, increase order by one")
    order <- order + 1
  } 
  if ( length(theta) != length(phi)) stop("getsphericalharmonicseven: need same length of theta and phi")
  kseq <- seq( 0, order, 2)
  n <- length( phi)
  values <- matrix( 0, (order+1)*(order+2)/2, n)
  #  if ( require( gsl, quietly = TRUE, warn.conflicts = FALSE)) {
  for (k in kseq) {
    mseq <- seq( -k, k, 1)
    for (m in mseq) {
      ind <- (k^2+k+2)/2+m
      z <- legendre_sphPlm( k, abs(m), cos(theta))## gsl
      if (m < 0) {
        z <- sqrt(2)*z*cos(m*phi)
      } 
      if (m > 0) {
        z <- sqrt(2)*z*sin(m*phi)
      }
      values[ ind, ] <- z
    }
  }
  #  } else {
  #    warning("gsl package not available \n returning zeros instead of estimates")
  #  }
  values
}


