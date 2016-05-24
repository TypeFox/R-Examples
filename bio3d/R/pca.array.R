"pca.array" <- function(x, use.svd=TRUE, ...) {
  if(!is.array(x))
    stop("provide an array of matrices")
  
  ## Log the call
  cl <- match.call()
  
  ## Construct the input matrix for PCA
  x <- t(apply(x, 3, function(y) y[upper.tri(y)]))
  
  dx <- dim(x)
  n <- dx[1]; p <- dx[2]
  mean <- colMeans(x)
  
  if(!use.svd) {
    ## coverance matrix
    S    <- var(x)          
    
    ## eigenvectors ("U") & eigenvalues ("L"): [ U'SU=L ]
    prj  <- eigen(S, symmetric = TRUE)
    L <- prj$values
    U <- prj$vectors
  }
  else {
    ## S = Q'Q, Q = UDV'
    Q <- t(t(x) - mean) / sqrt(n-1)
    prj <- svd(Q)
    L <- prj$d^2
    U <- prj$v
  }
  
  ## fix negative eigenvalues
  ## (these are very small numbers and should be zero)
  L[L<0]<-0
  sdev <- sqrt(L)
  
  ## scores of "xyz" on the pc's [ z=U'[x-x.mean] ]
  z <- sweep(x,2,mean) %*% (U)
  class(U)="pca.loadings"
  
  out <- list(L=L, U=U, z=z,
              sdev=sdev, mean=mean, call=cl)
  
  class(out)="pca"
  return(out)
}
