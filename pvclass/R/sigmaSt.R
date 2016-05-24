sigmaSt <-
function(X, Y, L, dimension, n, mu = NULL){
  if(any(is.null(mu))){
    # Compute mu
    mu <- matrix(0, L, dimension)
    for(b in seq_len(L)) {
      mu[b, ] = colMeans(X[Y == b, ])
    }
  }
  
  #Compute sigma
  Xc <- X - mu[Y, ]
  sigma <- crossprod(Xc) / (n - L)
}

