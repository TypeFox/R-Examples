dpostb <-
function(X,Z,Y,betas,gammas,bpri,Bpri) {
  
 
  
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
  
  if(!is.vector(bpri)){
    stop("The mean of the parameters must be a vector")
  }
  
  if(ncol(Bpri)!=nrow(Bpri)){
    stop("The covariance matriz is not square")
  }
  
#  if(ncol(X)!=nrow(betas)| ncol(Z)!=nrow(gammas)|ncol(X)!=length(b_pri)|ncol(X)!=ncol(B_pri)){
 #   stop("The initial values are not conformable to the data")
  #}
  
  eta <- X%*%betas
  mu <- exp(eta)/(1+exp(eta))
  
  tau <- Z%*%gammas
  phi <- exp(tau)
  
  L <- prod(dbeta(Y, phi*mu,phi*(1-mu))) # verosimilitud
  P <- dmvnorm(t(betas),bpri,Bpri) # priori
  value=L*P
  value
}
