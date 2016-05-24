#####################################################################################
# Update functions for the parameters which for an explicit solution does not exist #
# Solutions are calculated inside EM using root-finding algorithms                  #
#####################################################################################

# Update function for p.loggeom and theta
# ---------------------------------------
loggeom.update <- function(par, j, M, eta){
  p     <- plogis(par[1])
  theta <- plogis(par[2])
  a     <- eta[,j]
  b     <- log(dloggeom(1:M, p, theta))
  x     <- -sum(a[(a > 0)] * b[(a > 0)])
  return(x)
  }
