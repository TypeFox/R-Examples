qvolkov<-function(p, theta, m, J, lower.tail = TRUE, log.p = FALSE){
  if (length(theta) > 1 | length(m) > 1 |length(J) > 1) stop("vectorization of parameters is not implemented")
  if (!all(is.finite(c(J, theta, m))) | theta <= 0 | J <= 0 | m < 0 | m > 1)
	  return(rep(NaN, length(p)))
  if (log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  Y <- 1:J
  X <- pvolkov(Y, theta=theta, m=m, J=J)
  approx(X,Y,xout=p,method="constant",f=0,yleft=1, yright=J)$y
}
