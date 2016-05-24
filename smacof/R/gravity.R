gravity <- function(X, lambda = 1)
{
  Y <- X
  Y[Y > 0] <- 1
  co.maty <- t(Y)%*%Y       #co-occurences 
  grav.disty <- sqrt(diag(co.maty)%*%t(diag(co.maty))/co.maty) # gravity distance
  rownames(grav.disty) <- colnames(grav.disty)
  grav.disty1 <- grav.disty
  vec.disty <- grav.disty[lower.tri(grav.disty)][which(grav.disty[lower.tri(grav.disty)] != Inf)]
  W <- matrix(1, ncol = ncol(grav.disty), nrow = nrow(grav.disty))
  W[grav.disty1 == Inf] <- 0               #blank out Inf distances in SMACOF
  grav.disty[grav.disty1 == Inf] <- 1000   #replace Inf by any number (no matter which, they will be blanked out anyway)
  
  result <- list(gravdiss = grav.disty^lambda, weightmat = as.dist(W), co.occ = co.maty)
  result         
}