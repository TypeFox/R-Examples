
gillespied <- function (N, T=100, dt=1, ...)
{
  tt = 0
  n = T%/%dt
  x = N$M
  S = t(N$Post-N$Pre)
  u = nrow(S)
  v = ncol(S)
  xmat = matrix(ncol=u,nrow=n)
  i = 1
  target = 0
  repeat {
    h = N$h(x, tt, ...)
    h0 = sum(h)
    if (h0 < 1e-10)
      tt = 1e99
    else
      tt = tt+rexp(1,h0)
    while (tt >= target) {
      xmat[i,] = x
      i = i+1
      target = target+dt
      if (i > n)
        return(ts(xmat,start=0,deltat=dt))
    }
    j = sample(v,1,prob=h)
    x = x+S[,j]
  }
}


# eof

