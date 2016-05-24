gammaproposal <-
function(Y,X,Z,betas,gammas,gpri,Gpri){
  
  if(is.null(X)|is.null(Z)){
    stop("There is no data")
  }
  
  if(nrow(X)!=nrow(Z)){
    stop("The mean and precision data have a different size")
  }
  
  #if(!is.matrix(betas)|!is.matrix(gammas)){
   # stop("The parameters must be a matrix")
  #}
  
  if(!is.vector(gpri)){
    stop("the initial parameters for gamma must be a vector")
  }
  
  if(ncol(Gpri)!=nrow(Gpri)){
    stop("The initial covariance matrix for gamma is not square")
  }
  
 # if(ncol(X)!=nrow(betas)| ncol(Z)!=nrow(gammas)|ncol(Z)!=length(g_pri)|ncol(Z)!=ncol(G_pri)){
  #  stop("The initial values are not conformable to the data")
  #}  
  eta <- X%*%betas
  mu <- exp(eta)/(1+exp(eta)) 
  
  tau <- Z%*%gammas
  phi <- exp(tau)
  
  sigma <- (1-mu)/(mu*(1+phi))
  
  Y.tilde <- tau + Y/mu -1 
  
  G.pos <- (solve(solve(Gpri)+ t(Z)%*%solve(diag(as.vector(sigma)))%*%Z))
  g.pos <- G.pos%*%(solve(Gpri)%*%gpri + t(Z)%*%solve(diag(as.vector(sigma)))%*%Y.tilde)
  
  value=rmvnorm(1,g.pos,G.pos)
  value
}
