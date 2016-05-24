"HPloglik" <- function(X, theta, M, F, s, log=TRUE){
 val <- .Call("HPloglik", deltat(X), as.numeric(X), theta, M[[1]], M[[2]],
     M[[3]], M[[4]], M[[5]], M[[6]], M[[7]], F, s, .GlobalEnv, PACKAGE="sde")  
 if(log)
  return(val)

 exp(val)
}


"EULERloglik" <- function(X, theta, d, s, log=TRUE){
 val <-  .Call("EULERloglik", deltat(X), as.numeric(X), theta, d, s,
    .GlobalEnv, PACKAGE="sde") 
 if(log)
  return(val)

 exp(val)
}

SIMloglik <- function(X, theta, d, s,  M=10000, N=2, log=TRUE){
 val <- .Call("SIMloglik", as.numeric(X), deltat(X), d, s, theta, as.integer(N), 
  as.integer(M), .GlobalEnv, PACKAGE = "sde")
 if(log)
  return(val)
 exp(val) 
}

