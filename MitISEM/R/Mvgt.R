# Sampling and density evaluation for the mixture of p-variate 
# student t densities
#
# inputs:
#    N     : [double] number of draws
#    mit   : [list] parameters of the mixture t density (see 'is.Mit')
#    theta : [matrix Nxp] matrix of draws
#    log   : [logical] 'dmvgt' returns log-density values if 'log=TRUE'
# outputs:
#    theta : [matrix Nxp] draws from the density
#    dens  : [vector size N] (log) density at parameter values
#
# author : Nalan Basturk
# date   : 20120912

# random sampling from mixture of t densities
rmvgt <- function(N = 1, mit = list()){
  # test if mit is well defined
  if(!isMit(mit)){
    warning("'mit' not well defined in 'rmvgt'; set to univariate std. student t")
    mit <- list(p = 1, mu = as.matrix(0), Sigma = as.matrix(1),df = 1)
  }
  H     <- length(mit$p)        # number of components
  k     <- ncol(mit$mu)         # dimension of t dist.
  # sample memberships
  memb  <- sample(1:H, N, prob = mit$p, replace = TRUE)
  theta <- matrix(nrow=N, ncol=k)
  for(h in 1:H){
    ind_h  <- (memb==h)
    n_h    <- sum(ind_h)
    if(n_h > 0){
      mu_h   <- mit$mu[h,]
      draw_h <- rmvt(n_h, matrix(mit$Sigma[h,],k,k), mit$df[h])
      draw_h  = draw_h + matrix(mu_h, n_h, k, byrow=TRUE)
      theta[ind_h,] = draw_h
    }
  }    
  return(theta)
}

# density of a multivariate t density
dmvgt <- function(theta, mit = list(), log = TRUE){
  if (missing(theta)) 
    stop("'theta' is missing in 'dMit'")  
  if(is.vector(theta)) 
    theta <- matrix(theta,nrow=1)
  # test if mit is well defined
  if(!isMit(mit)){ 
    mit <- list(p = 1, mu = as.matrix(0), Sigma = as.matrix(1),df = 1)
    warning("'mit' not well defined in 'dmvgt'; set to univariate std. student t")
  }
  H     <- length(mit$p)        # number of components
  k     <- ncol(mit$mu)         # dimension of t dist.
  theta  = as.matrix(theta,nrow=k)
  N     <- nrow(theta)   
  dcoms <- matrix(nrow=N,ncol=H)
  for(h in 1:H)
    dcoms[,h] = dmvt(theta,mit$mu[h,],matrix(mit$Sigma[h,],k,k),mit$df[h],log=TRUE)
  tmp  <- log(matrix(mit$p,N,H,byrow=TRUE)) + dcoms
  dens <- rowSums(exp(tmp))  
  if (log) 
    dens <- log(dens)
  if (!log)
    dens <- dens
  return(dens)
}
