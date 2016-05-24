qgs<-function(p, k, S, lower.tail = TRUE, log.p = FALSE){
  if (length(S) > 1) stop("vectorization of S is not implemented")
  if (S <= 0 | !is.wholenumber(S) ) return(rep(NaN, length(p)))
  if (k<0 | k>1) return(rep(NaN, length(p)))
  if (log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  ## Ugly: just to make qgs(1, ...) = S
  p[p==1] <- 1+1e-12
  Y <- 1:S
  X <- pgs(Y, k=k, S=S)
  approx(x=X, y=Y, xout=p, method="constant", rule=2)$y
}
