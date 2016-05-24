modwt.backward <- function(W, V, filter, j){

  # error checking
  if(class(W) != "numeric")
    stop("Invalid argument: 'W' must be of class 'numeric' or 'matrix'.")
  if(class(V) != "numeric")
    stop("Invalid argument: 'V' must be of class 'numeric' or 'matrix'.")
  if(class(filter) != "wt.filter")
    stop("Invalid argument: 'filter' must be of class 'wt.filter'.")
  if(filter@transform == "dwt")
    stop("Invalid argument: 'transform' element of agrument 'filter' must be equal to 'modwt'.")
  if(is.na(match(class(j), c("numeric", "integer"))))
    stop("Invalid argument: 'j' must be of class 'numeric' or 'integer'.")
  if(length(j) != 1)
    stop("Invalid argument: 'j' must have length 1.")

  N <- length(V)
  h <- filter@h
  g <- filter@g
  L <- length(h)
  Vj <- rep(NA, length=N)
  for(t in 0:(N-1)){
    k <- t
    Vjt <- h[1]*W[k+1] + g[1]*V[k+1]
    for(n in 1:(L-1)){
      k <- k + 2^(j-1)
      if(k >= N) k <- k - floor(k/N)*N
      Vjt <- Vjt + h[n+1]*W[k+1] + g[n+1]*V[k+1]
    }
    Vj[t+1] <- Vjt
  }
  return(Vj)
}
