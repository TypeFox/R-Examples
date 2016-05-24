# Note:  normmixrm.sim is here for backwards compatibility
rmvnormmix <- normmixrm.sim <- function(n,lambda=1,mu=0,sigma=1) {
  m <- length(lambda) # nb of components
  mu <- matrix(mu, nrow=m)
  sigma <- matrix(sigma, nrow=m)
  if ((r <- NCOL(mu)) != NCOL(sigma)) {
    stop("mu and sigma must have the same number of columns", call.=FALSE)
  }
  z <- sample(m,n,replace=TRUE,prob=lambda) # component 
  matrix(rnorm(n*r,mean=as.vector(mu[z,]),sd=as.vector(sigma[z,])),n,r)  
}

