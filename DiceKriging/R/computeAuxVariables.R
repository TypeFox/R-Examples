compute.beta.hat <- function(x, M){
  l <- lm(x ~ M - 1)
  return(as.numeric(l$coef))
}

compute.z <- function(x, M, beta=NULL){
  if (!is.null(beta)) {
    z <- x - M %*% beta
  } else {   # in that case, z depends on the MLE estimate of beta (which is useless to compute here)
    Q <- qr.Q(qr(M))
    H <- Q %*% t(Q)
    z <- x - H %*% x
  }
  return(as.numeric(z))
}

compute.sigma2.hat <- function(z){
  drop(crossprod(z)/length(z))
}


computeAuxVariables <- function(model) {
	  
	aux <- covMatrix(model@covariance, X=model@X, noise.var=model@noise.var)
	C <- aux[[1]]
  T <- chol(C)
	
  x <- backsolve(t(T), model@y, upper.tri = FALSE)
  M <- backsolve(t(T), model@F, upper.tri = FALSE)
  model@T <- T
	model@M <- M
	
  if (length(model@trend.coef)>0){
    z <- backsolve(t(T), model@y-model@F%*%as.matrix(model@trend.coef), upper.tri=FALSE)   
    model@z <- as.numeric(z)
  } else {
    model@z <- numeric(0)
  }
	
  return(model)
} 
