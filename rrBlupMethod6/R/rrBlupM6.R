rrBlupM6 <- function(y,
                     X = matrix(1,nrow=n,ncol=1),
                     Z,
                     sig2e,
                     chunks = as.integer(1)){
  
  p <- ncol(Z) ## no. of markers
  n <- length(y) ## no. of observations
  

  Zt <- t(Z)
  Xt <- t(X)

  tolerance <- .Machine$double.eps
  bounds <- c(1e-09, 1e+09)

  

  if(identical(as.integer(chunks),as.integer(1))){gamma <- Z %*% t(Z)}
  else{
    ## Compute gamma from chunks
    chunksize <- floor(p/chunks)
    check <- chunks*chunksize-p
    gamma <- matrix(0, nrow(Z), nrow(Z))
    gamma <- as.matrix(gamma)
    
    if (check<0) { chunks=chunks+1}
    
    for (i in 1:chunks) {
      from <- (i-1)*chunksize+1
      to <- i*chunksize
      if (to>p) {to=p}
      Zti <- Zt[from:to, ,drop = FALSE]
      gamma <- gamma + t(Zti)%*%Zti
    }
  }
  

  
  ## Get eta and lambda via spectral decomposition of SHbS
  S <- diag(n) - (X %*% solve(t(X) %*% X) %*% t(X))     
  SHbS <- S %*% (gamma + diag(n)) %*% S
  SVS <- eigen(SHbS)
  lambda_asMatrix <- SVS$values - 1
  lambda_asMatrix <- as.matrix(lambda_asMatrix)
  vec <- 1:n

  ## get non-zero eigenvalues
  for( i in vec ){
    if(lambda_asMatrix[i:i,]>1e-15){
      lambda <- lambda_asMatrix[1:i,]
      r <- i
    }
  }
  
  Q <- SVS$vectors 
  eta <- crossprod(Q, y)
  eta.sq <- as.matrix(eta^2)
  eta.sq <- eta.sq [1:r,]

  f.REML <- function(sig2u,
                     sig2e,
                     lambda,
                     eta.sq) {
    sum(eta.sq/(sig2u * lambda + sig2e)) +
      sum(log(sig2u * lambda + sig2e))   
  }

  soln <- optimize(f.REML,
                   interval = bounds,
                   sig2e,
                   lambda,
                   eta.sq,
                   tol = tolerance)
  
  sig2u <- soln$minimum

  ## Compute inverse of V and estimate of fixed parameter
  Diag_Se <- diag(n) * sig2e
  Hinv <- solve(gamma * sig2u + Diag_Se)

  ## BLUE of fixed effects
  betahat <- solve(Xt %*% Hinv %*% X, Xt %*% Hinv %*% y)

  ## BLUP of u
  uhat <- sig2u * Zt %*% Hinv %*% (y - X %*% betahat)

  return(list(uhat = as.vector(uhat),
              betahat = as.vector(betahat),
              sig2u = sig2u)
         )
}



