dcSim <- function(x0, x, t, d, s, theta, M=10000, N=10, log=FALSE){
 val <- .Call("dcSim", x0, x, t, d, s, theta, as.integer(N), 
  as.integer(M), .GlobalEnv, PACKAGE = "sde")
 if(log)
  return(log(val))
 val 
}

