qmand <- function (p, N, s, v, lower.tail = TRUE, log.p = FALSE) {
  if (length(N) > 1) stop("vectorization of N is not implemented")
  if (N < 1 |!is.wholenumber(N)) return(rep(NaN, length(p)))
  if (log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  ## Ugly: just to make qmand(1, ...) = S
  p[p==1] <- 1+1e-12
  Y <- 1:N
  X <- pmand(Y, N=N, s=s, v=v)
  approx(x=X, y=Y, xout=p, method="constant", rule=2)$y
}
