####################################################
## Random generation from a contaminated Gaussian ##
####################################################

rCN <- function(n, mu = rep(0,p), Sigma, alpha = 0.99, eta = 1.01){
  
  if(missing(Sigma))
    stop("Sigma is missing")
  if(alpha<0 | alpha>1)
    stop("alpha must be in (0,1)")
  if(eta<1)
    stop("eta must be greater than 1")
  
  p <- if(is.matrix(Sigma)) 
    ncol(Sigma)
  else 1
  
  X    <- array(0,c(n,p),dimnames=list(1:n,paste("X.",1:p,sep="")))
  good <- rbinom(n=n,size=1,prob=alpha)
  for(i in 1:n){
    if(good[i]==1)
      X[i,] <- rmnorm(n = 1, mean = mu, varcov = Sigma)
    else
      X[i,] <- rmnorm(n = 1, mean = mu, varcov = eta*Sigma)
  }
  return(X)  
  
}
