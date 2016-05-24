working.comp <- function(fit){
   
  G <- -(fit$gradient - fit$S.h2)
  H <- fit$hessian - fit$S.h
  
  W.eig <- eigen(H, symmetric=TRUE)
  
  if(min(W.eig$values) < sqrt(.Machine$double.eps)) { pep <- which(W.eig$values < sqrt(.Machine$double.eps)); W.eig$values[pep] <- 0.0000001 }  
  
  srev    <- sqrt(W.eig$val)
  
  if(length(srev) == 1) {
    c.W     <- W.eig$vec%*%tcrossprod(srev, W.eig$vec) 
    W.invsr <- W.eig$vec%*%tcrossprod(1/srev, W.eig$vec)
  } else {
    c.W     <- W.eig$vec%*%tcrossprod(diag(srev)  ,W.eig$vec) 
    W.invsr <- W.eig$vec%*%tcrossprod(diag(1/srev),W.eig$vec)
  }
  
  X <- c.W 
  Z <- W.invsr%*%G + X%*%fit$argument 
                 
 list( X = X , Z = Z )

}


























