"pwm" <-
function(x,nmom=5,sort=TRUE) {
  if(sort) x <- sort(x)
  n <- length(x)

  betas <- vector(mode="numeric")
  for(r in seq(0,nmom-1)) {
    i <- r+1
    sum <- 0
    for(j in seq(1,n)) {
      sum <- sum + choose(j-1,r)*x[j]
    }
    betas[i] <- sum/(n*choose(n-1,r))
  }
  z <- list(betas=betas,source="pwm")
  return(z)
}


# D <- rnorm(10)
# pwm(D)
# lmom2pwm(lmorph(lmoms(D)))
# lmom2pwm(lmom.ub(D))
