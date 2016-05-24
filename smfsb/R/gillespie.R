
gillespie <- function(N, n, ...)
{
  tt = 0
  x = N$M
  S = t(N$Post-N$Pre)
  u = nrow(S)
  v = ncol(S)
  tvec = vector("numeric",n)
  xmat = matrix(ncol=u,nrow=n+1)
  xmat[1,] = x
  for (i in 1:n) {
    h = N$h(x,tt, ...)
    tt = tt+rexp(1,sum(h))
    j = sample(v,1,prob=h)
    x = x+S[,j]
    tvec[i] = tt
    xmat[i+1,] = x
  }
  return(list(t=tvec, x=xmat))
}

discretise <- function(out, dt=1, start=0)
{
  events = length(out$t)
  end = out$t[events]
  len = (end-start)%/%dt+1
  x = matrix(nrow=len,ncol=ncol(out$x))
  target = 0
  j = 1
  for (i in 1:events) {
    while (out$t[i] >= target) {
      x[j,] = out$x[i,]
      j = j+1
      target = target+dt
    }
  }
  ts(x, start=0, deltat=dt)
}


# eof

