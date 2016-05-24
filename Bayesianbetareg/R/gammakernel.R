gammakernel <-
function(X,Z,Y,gammas.n,betas.v,gammas.v,gpri,Gpri){
 
  
  Y=as.matrix(Y)
  
  if(is.null(X)|is.null(Z)|is.null(Y)){
    stop("There is no data")
  }
  
  if(nrow(X)!=nrow(Z)|nrow(X)!=nrow(Y)|nrow(Z)!=nrow(Y)){
    stop("The variables have a diferent size")
  }
  
  #if(!is.matrix(gammas.n)|!is.matrix(gammas.v)|!is.matrix(betas.v)){
   # stop("The parameters must be a matrix")
  #}
  
  if(!is.vector(gpri)){
    stop("The mean of the parameters must be a vector")
  }
  
  if(ncol(Gpri)!=nrow(Gpri)){
    stop("The covariance matriz is not square")
  }
  
 # if(ncol(X)!=nrow(betas.v)|ncol(Z)!=nrow(gammas.v)| ncol(Z)!=nrow(gammas.n)|ncol(Z)!=length(g_pri)|ncol(Z)!=ncol(G_pri)){
  #  stop("The initial values are not conformable to the data")
  #}   
  eta <- X%*%betas.v
  mu <- exp(eta)/(1+exp(eta)) 
  
  tau <- Z%*%gammas.v
  phi <- exp(tau)
  
  sigma <- (1-mu)/(mu*(1+phi))
  
  Y.tilde <- tau + Y/mu -1 
  
  G.pos <- solve(solve(Gpri)+ t(Z)%*%solve(diag(as.vector(1/sigma)))%*%Z)
  g.pos <- G.pos%*%(solve(Gpri)%*%gpri + t(Z)%*%solve(diag(as.vector(1/sigma)))%*%Y.tilde)
  
  value=dmvnorm(t(gammas.n),g.pos,G.pos)
  value
}
