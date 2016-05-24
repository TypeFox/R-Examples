muproposal <-
function(Y,X,Z,betas,gammas,bpri,Bpri){
  
  
  
  
  if(is.null(X)|is.null(Z)){
    stop("There is no data")
  }
  
  if(nrow(X)!=nrow(Z)){
    stop("The mean and precision data have a different size")
  }
  
  #if(!is.matrix(betas)|!is.matrix(gammas)){
   # stop("The parameters must be a matrix")
  #}
  
  if(!is.vector(bpri)){
    stop("the initial parameters for beta must be a vector")
  }
  
  if(ncol(Bpri)!=nrow(Bpri)){
    stop("The initial covariance matrix for beta is not square")
  }
  
  #if(ncol(X)!=nrow(betas)| ncol(Z)!=nrow(gammas)|ncol(X)!=length(b_pri)|ncol(X)!=ncol(B_pri)){
   # stop("The initial values ??are not conformable to the data")
  #}
  
  eta <- X%*%betas
  mu  <- exp(eta)/(1+exp(eta)) 
  
  tau <- Z%*%gammas
  phi <- exp(tau)
  
  Y.mono <- eta + (Y-mu)/(mu*(1-mu))
  sigmab <- 1/((mu*(1-mu))*(1+phi))
  
  B.pos <- (solve(solve(Bpri)+ t(X)%*%solve(diag(as.vector(sigmab)))%*%X))
  b.pos <- B.pos%*%(solve(Bpri)%*%bpri + t(X)%*%solve(diag(as.vector(sigmab)))%*%Y.mono)
  
  value=rmvnorm(1,b.pos,B.pos)
  value
}
