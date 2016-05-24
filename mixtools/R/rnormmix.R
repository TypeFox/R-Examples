# Simulate from a normal mixture.  New version simply calls rmvnormmix.

# normmix.sim is here for backwards compatibility
rnormmix <- normmix.sim <- function(n,lambda=1,mu=0,sigma=1) {
  if (NCOL(mu)>1 || NCOL(sigma)>1)
    stop ("Use the rmvnormmix function instead.")
  as.vector(rmvnormmix(n,lambda,mu,sigma))
}

