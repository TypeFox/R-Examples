qbs<-function(p, N, S, lower.tail = TRUE, log.p = FALSE){
  if (length(N) > 1 | length(S) > 1) stop("vectorization of parameters is not implemented")
  if (N <= 0 | S <= 0 | !is.wholenumber(N) | !is.wholenumber(S))  return(rep(NaN, length(p)))
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  ## Ugly: just to make qbs(1, ...) = N
  p[p==1] <- 1+1e-12
  Y <- 1:N
  X <- pbs(Y, N=N, S=S)
  approx(x=X, y=Y, xout=p, method="constant", rule=2)$y
}
