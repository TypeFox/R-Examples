##
## Dirichlet (Multivariate Beta)
##
##  From MCMCpack 0.7-2 by Andrew Martin and Kevin Quinn
#
#  
# note: this code is taken verbatim from the R-package
# "Greg's Miscellaneous Functions" maintained by
# Gregory R. Warnes <Gregory_R_Warnes@groton.pfizer.com>
#
# Kevin Rompala 5/6/2003
rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE)
  sm <- x%*%rep(1,l)
  return(x/as.vector(sm))
}
