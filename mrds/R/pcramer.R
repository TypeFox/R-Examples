pcramer <-
function (q, eps = 1e-05) 
{
#
#  Function taken from coda package
#
  log.eps <- log(eps)
  y <- matrix(0, nrow = 4, ncol = length(q))
  for (k in 0:3) {
    z <- gamma(k+0.5)*sqrt(4*k+1)/(gamma(k+1)*pi^(3/2)*sqrt(q))
    u <- (4*k+1)^2/(16*q)
    y[k+1, ] <- ifelse(u >-log.eps, 0, z * exp(-u)*besselK(x=u,nu=1/4))
  }
  return(apply(y, 2, sum))
}
