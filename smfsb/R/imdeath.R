
imdeath <- function(n=20, x0=0, lambda=1, mu=0.1)
{
  xvec = vector("numeric",n+1)
  tvec = vector("numeric",n)
  t = 0
  x = x0
  xvec[1] <- x
  for (i in 1:n) {
    t = t+rexp(1,lambda+x*mu)
    if ( runif(1,0,1) < lambda/(lambda+x*mu) )
      x <- x+1
    else
      x <- x-1
    xvec[i+1] <- x
    tvec[i] <- t
  }
  stepfun(tvec, xvec)
}

# eof

