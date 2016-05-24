hr05AdjustedDF <- 
#
# Calculate the adjusted degrees of freedom parameter
# for the F distribution given in Hardin and Rocke (2005)
# for testing Mahalanobis distances calculated with the 
# MCD
# 
# Christopher G. Green
# 2011
#
function( n.obs, p.dim, 
  mcd.alpha = max.bdp.mcd.alpha(n.obs, p.dim), 
  m.asy = ch99AsymptoticDF( n.obs, p.dim, mcd.alpha )$m.hat.asy, 
  method=c("HR05","GM14") ) 
{

  method <- match.arg(method)
  retval <- numeric(0)

  retval <- if ( method == "HR05" ) {
    if ( mcd.alpha == max.bdp.mcd.alpha(n.obs, p.dim) ) {
      # original equation from Hardin and Rocke 2005
      hr05.predict.050.hr05(m.asy,p.dim,n.obs)
    } else {
      stop("HR05 unsupported for alpha other than maximum breakdown case.")
    }
  } else if ( method == "GM14" ) {
    # use fitted models to adjust asymptotic degrees of freedom to 
    # simulated values for small samples
    predictfunc( m.asy, p.dim, n.obs, mcd.alpha )
  }

  retval

}

# internal functions

hr05.predict.050.hr05 <- function(m.asy, p, n)  { 
  # original equation from Hardin and Rocke 2005
  m.asy * exp( 0.725  - 0.00663*p - 0.0780*log(n)) 
}

predictfunc <- function(m.asy, p, n, alpha) {
  # use fitted models to adjust asymptotic degrees of freedom to 
  # simulated values for small samples
  #z1 <- 13.1072594 - 15.0079984 * alpha + 0.1347435 * p
  #z2 <- n^(0.5119686 + 0.2235726 * alpha)
  # 2014-03-01 update
  #z1 <- 13.263331 - 15.092671 * alpha + 0.126725 * p
  #z2 <- n^(0.567457 + 0.149904 * alpha)
  # 2014-03-29 update
  z1 <- 12.745653 - 14.545559 * alpha + 0.127400 * p
  z2 <- n^(0.559217 + 0.149040 * alpha)
  m.asy * exp( z1/z2 )
}

