graphicalVARsim <- function(
  nTime, # Number of time points
  beta,
  kappa,
  init = 0,
  intercepts = 0,
  warmup = 100){
  
  stopifnot(!missing(beta))
  stopifnot(!missing(kappa))  
  
  Nvar <- ncol(kappa)
  init <- rep(init, length = Nvar)
  intercepts <- rep(intercepts, length = Nvar)
  
  totTime <- nTime + warmup
  
  Data <- matrix(NA, totTime, Nvar)
  Data[1,] <- init
  
  Sigma <- solve(kappa)

  for (t in 2:totTime){
    Data[t,] <- t(intercepts + beta %*% Data[t-1,]) + rmvnorm(1, rep(0,Nvar), Sigma)
  }
  
  return(Data[-seq_len(warmup), ,drop=FALSE])
}