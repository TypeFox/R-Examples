dpostg <-
function(X,Z,Y,betas,gammas,gpri,Gpri) {
  
  
  
  Y=as.matrix(Y)
  
  if(is.null(X)|is.null(Z)|is.null(Y)){
    stop("There is no data")
  }
  
  if(max(Y)> 1 | min(Y) < 0){
    stop("The data is not in the (0,1) interval")
  }
  
  if(nrow(X)!=nrow(Z)|nrow(X)!=nrow(Y)|nrow(Z)!=nrow(Y)){
    stop("The variables have a diferent size")
  }
  
  #if(!is.matrix(betas)|!is.matrix(gammas)){
   # stop("The parameters must be a matrix")
  #}
  
  if(!is.vector(gpri)){
    stop("The mean of the parameters must be a vector")
  }
  
  if(ncol(Gpri)!=nrow(Gpri)){
    stop("The covariance matriz is not square")
  }
  
#  if(ncol(X)!=nrow(betas)| ncol(Z)!=nrow(gammas)|ncol(Z)!=length(g_pri)|ncol(Z)!=ncol(G_pri)){
 #   stop("The initial values are not conformable to the data")
  #} 
  eta <- X%*%betas
  mu <- exp(eta)/(1+exp(eta))
  
  tau <- Z%*%gammas
  phi <- exp(tau)
  
  L <- prod(dbeta(Y, phi*mu,phi*(1-mu))) # verosimilitud
  P <- dmvnorm(t(gammas),gpri,Gpri) # priori
  value=L*P
  value
}
