expsnps <-
function(x, m.rate, c.rate, tau) {
  veclen <- length(x)
  sumterm <- numeric(veclen)
  for (k in 1:veclen) {
    for (i in 0:x[k]) {
      sumterm[k] <- sumterm[k] + dgeom(i, c.rate/(2*m.rate+c.rate))*dpois(x[k]-i, m.rate*tau)
    }
  }
  return(sumterm)
}
