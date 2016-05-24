
rdiff <- function(afun, bfun, x0 = 0, t = 50, dt = 0.01, ...)
{
  n <- t/dt
  xvec <- vector("numeric", n)
  x <- x0
  sdt <- sqrt(dt)
  for (i in 1:n) {
    t <- i*dt
    x <- x + afun(x,...)*dt + 
      bfun(x,...)*rnorm(1,0,sdt)
    xvec[i] <- x
  }
  ts(xvec, deltat = dt)
}


# eof




