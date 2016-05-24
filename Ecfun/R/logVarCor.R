logVarCor <- function(x, corr, ...){
  # logVarCor represents a non-negative definite 
  # matrix as log(diag(x)) with corr attribute 
  #  
  # if matrix, return log(diag(x)) with attr(, 'corr')
  dimx <- dim(x)
  if(length(dimx)>1){
    if(length(dimx)>2){
      stop('x is not a matrix; dim = ', 
           paste(dimx, collapse=', '))
    }
    logx <- log(diag(x))
    names(logx) <- rownames(x)
    Corr <- cov2cor(x)
    if(dimx[1]<2){
      attr(logx, 'corr') <- numeric(0)
    } else {
      attr(logx, 'corr') <- Corr[lower.tri(Corr)]
    }
    return(logx)
  }
  # otherwise, construct matrix from lower.tri(corr)
  # with diag = exp(x)
  k <- length(x)
  if(missing(corr)){
    corr <- attr(x, 'corr')
    if(length(corr)<1){
      corr <- rep(0, choose(k, 2))
    }
  }  
  s <- exp(x/2)
  S2 <- outer(s, s)
  if(k<2){
    X <- S2 
    return(X)
  }
  Cor1 <- matrix(0, k, k)
  Cor1[lower.tri(Cor1)] <- corr
  Cor <- (diag(k)+Cor1+t(Cor1))
  X <- (S2*Cor)
  X
}