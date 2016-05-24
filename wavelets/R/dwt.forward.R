dwt.forward <- function(V, filter){
  
  # error checking
  if(class(V) != "numeric")
    stop("Invalid argument: 'V' must be of class 'numeric' or 'matrix'.")
  if(class(filter) != "wt.filter")
    stop("Invalid argument: 'filter' must be of class 'wt.filter'.")
  if(filter@transform == "modwt")
    stop("Invalid argument: 'transform' element of agrument 'filter' must be equal to 'dwt'.")

  h <- filter@h
  g <- filter@g
  L <- filter@L
  M <- length(V)
  Wj <- rep(NA, length=(M/2))
  Vj <- rep(NA, length=(M/2))
  for(t in 0:(M/2 - 1)){
    u <- 2*t + 1
    Wjt <- h[1]*V[u+1]
    Vjt <- g[1]*V[u+1]
    for(n in 1:(L-1)){
      u <- u - 1
      if(u < 0) u <- M - 1
      Wjt <- Wjt + h[n+1]*V[u+1]
      Vjt <- Vjt + g[n+1]*V[u+1]
    }
    Wj[t+1] <- Wjt
    Vj[t+1] <- Vjt
  }
  results <- list(W = Wj, V = Vj)
  return(results)
}
