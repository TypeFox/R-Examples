
rcfmc <- function(n,Q,pi0)
{
  xvec = vector("numeric",n+1)
  tvec = vector("numeric",n)
  r = length(pi0)
  x = sample(r,1,prob=pi0)
  t = 0
  xvec[1] = x
  for (i in 1:n) {
    t = t+rexp(1,-Q[x,x])
    weights = Q[x,]
    weights[x] = 0
    x = sample(r,1,prob=weights)
    xvec[i+1] = x
    tvec[i] = t
  }
  stepfun(tvec,xvec)
}


# eof
