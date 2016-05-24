BBridge <- function(x=0, y=0, t0=0, T=1, N=100){
  if(T<= t0) stop("wrong times")
  dt <- (T-t0)/N
  t <- seq(t0, T, length=N+1)
  X <- c(0,cumsum( rnorm(N)*sqrt(dt)))
  BB <- x + X - (t-t0)/(T-t0)*(X[N+1]-y+x)
  X <- ts(BB, start=t0,deltat=dt)
  return(invisible(X))
}

BM <- function(x=0, t0=0, T=1, N=100){
  if(T<= t0) stop("wrong times")
  dt <- (T-t0)/N
  t <- seq(t0,T, length=N+1)
  X <- ts(cumsum(c(x,rnorm(N)*sqrt(dt))),start=t0, deltat=dt)
  return(invisible(X))
}

GBM <- function(x=1, r=0, sigma=1, T=1, N=100){
   B <- BM(T=T,N=N)
   S <- x * exp((r-sigma^2/2)*time(B) + sigma* as.numeric(B))
   X <- ts(S, start=0,deltat=deltat(B))
   return(invisible(X))
}
