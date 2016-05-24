#######################################################
# Density functions of the sojourn time distributions #
#######################################################

# Density function of the logarithmic-geometric distribution
# ----------------------------------------------------------
dloggeom <- function(x, p, theta){
  # Auxiliary function
  f <- function(x, k){
    return(choose(k + x - 1, x) * (1 - p)^k * p^x * theta^k / k / (-log(1 - theta)))
    }
  lower    <- 1e-10
  k.max    <- 500
  finished <- as.logical(FALSE)
  C <- (f(aperm(matrix(x, byrow = TRUE, 1)) %*% matrix(rep(1, k.max), byrow = TRUE, 1),
   		  aperm(matrix(rep(1, length(x)), byrow = TRUE, 1)) %*% matrix(c(1:k.max), byrow = TRUE, 1)))
  C[is.na(C) | (C == Inf)] <- 0
  return(array(matrix(rep(1, k.max), byrow = TRUE, 1) %*% aperm(C)))
  }
