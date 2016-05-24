mukernel <-
function(X,Z,Y,betas.n,betas.v,gammas.v,bpri,Bpri){
 
  
  Y=as.matrix(Y)
  
  if(is.null(X)|is.null(Z)|is.null(Y)){
    stop("There is no data")
  }
  
  if(nrow(X)!=nrow(Z)|nrow(X)!=nrow(Y)|nrow(Z)!=nrow(Y)){
    stop("The variables have a diferent size")
  }
  
  #if(!is.matrix(betas.n)|!is.matrix(gammas.v)|!is.matrix(betas.v)){
   # stop("The parameters must be a matrix")
  #}
  
  if(!is.vector(bpri)){
    stop("The mean of the parameters must be a vector")
  }
  
  if(ncol(Bpri)!=nrow(Bpri)){
    stop("The covariance matriz is not square")
  }
  
  #if(ncol(X)!=nrow(betas.n)|ncol(X)!=nrow(betas.v)| ncol(Z)!=nrow(gammas.v)|ncol(X)!=length(b_pri)|ncol(X)!=ncol(B_pri)){
   # stop("The initial values are not conformable to the data")
  #}
  
  eta <- X%*%betas.v
  mu <- exp(eta)/(1+exp(eta)) 
  
  tau <- Z%*%gammas.v
  phi <- exp(tau)
  
  Y.mono <- eta + (Y-mu)/(mu*(1-mu))
  sigmab <- 1/((mu*(1-mu))*(1+phi))
  
  B.pos <- solve(solve(Bpri)+ t(X)%*%solve(diag(as.vector(sigmab)))%*%X)
  b.pos <- B.pos%*%(solve(Bpri)%*%bpri + t(X)%*%solve(diag(as.vector(sigmab)))%*%Y.mono)
  
  
  value=dmvnorm(t(betas.n),b.pos,B.pos)
  value
}
