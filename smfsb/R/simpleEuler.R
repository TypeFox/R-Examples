# function

simpleEuler <- function(t=50, dt=0.001, fun, ic, ...)
{
  p = length(ic)
  n = t/dt
  xmat = matrix(0,ncol=p,nrow=n)
  x = ic
  t = 0
  xmat[1,] = x
  for (i in 2:n) {
    t = t + dt
    x = x + fun(x,t,...) * dt
    xmat[i,] = x
  }
  ts(xmat, start=0, deltat=dt)
}


# eof

