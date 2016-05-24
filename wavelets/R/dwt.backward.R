dwt.backward <- function(W, V, filter){

  # error checking
  if(class(W) != "numeric")
    stop("Invalid argument: 'W' must be of class 'numeric' or 'matrix'.")
  if(class(V) != "numeric")
    stop("Invalid argument: 'V' must be of class 'numeric' or 'matrix'.")
  if(class(filter) != "wt.filter")
    stop("Invalid argument: 'filter' must be of class 'wt.filter'.")
  if(filter@transform == "modwt")
    stop("Invalid argument: 'transform' element of agrument 'filter' must be equal to 'dwt'.")

  M <- length(V)
  h <- filter@h
  g <- filter@g
  L <- length(h)
  Vj <- rep(NA, length=2*M)
  l <- -2
  m <- -1
  for(t in 0:(M-1)){
    l <- l+2
    m <- m+2
    u <- t
    i <- 1
    k <- 0
    Vj[l+1] <- h[i+1]*W[u+1] + g[i+1]*V[u+1]
    Vj[m+1] <- h[k+1]*W[u+1] + g[k+1]*V[u+1]
    if(L > 2){
      for(n in 1:((L/2)-1)){
        u <- u+1
        if(u >= M) u <- 0
        i <- i+2
        k <- k+2
        Vj[l+1] <- Vj[l+1] + h[i+1]*W[u+1] + g[i+1]*V[u+1]
        Vj[m+1] <- Vj[m+1] + h[k+1]*W[u+1] + g[k+1]*V[u+1]
      }
    }
  }
  return(Vj)
}
