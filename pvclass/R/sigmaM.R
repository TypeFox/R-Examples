sigmaM <-
function(X, Y, L, dimension, n, nvec, mu = NULL, sigma = NULL){
  if(is.null(mu)) {
    # Compute mu
    mu <- matrix(0, L, dimension)
    for(b in seq_len(L)) {
      mu[b, ] = colMeans(X[Y == b, ])
    }
  }
  
  M <- array(0, c(L, dimension, dimension))
  Xc <- X - mu[Y, ]
  for(b in seq_len(L)){
    M[b, , ] <- crossprod(Xc[Y == b, ])
  }

  if(is.null(sigma)) {
    # Compute sigma
    sigma <- colSums(M, dims = 1) / (n - L)
  }

  Mnvec <- nvec * M
  dimn <- dimension / n
  Id <- diag(1, dimension)
  
  # Fixed point iteration
  repeat{
    sigma.old <- sigma
    sigma.inv <- solve(sigma)
    sigma <- matrix(0, dimension, dimension)
    for(b in seq_len(L)){
      sigma <- sigma + Mnvec[b, , ] / sum(diag( sigma.inv %*% M[b, , ] ))
    }
    sigma <- dimn * sigma
    
    if(Matrix::norm(sigma.inv %*% sigma - Id, type = "F") < 10^-5) break
  }
  return(sigma)
}

