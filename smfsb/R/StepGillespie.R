# function

StepGillespie <- function(N)
{
  S = t(N$Post-N$Pre)
  v = ncol(S)
  return(
         function(x0, t0, deltat, ...)
         {
           t = t0
           x = x0
           termt = t0+deltat
           repeat {
             h = N$h(x,t,...)
             h0 = sum(h)
             if (h0 < 1e-10)
               t = 1e99
             else if (h0 > 1e6) {
               t = 1e99
               warning("Hazard too big - terminating simulation!")
             }
             else
               t = t+rexp(1,h0)
             if (t >= termt)
               return(x)
             j = sample(v,1,prob=h)
             x = x+S[,j]
           }
         }
         )		
}	

simTs <- function(x0, t0=0, tt=100, dt=0.1, stepFun, ...)
{
  n = (tt-t0) %/% dt + 1
  u = length(x0)
  names = names(x0)
  mat = matrix(nrow=n,ncol=u)
  x = x0
  t = t0
  mat[1,] = x
  for (i in 2:n) {
    t = t+dt
    x = stepFun(x,t,dt,...)
    mat[i,] = x
  }
  ts(mat, start=t0, deltat=dt, names=names)
}

simTimes <- function(x0, t0=0, times, stepFun, ...)
{
  n = length(times)
  u = length(x0)
  names = names(x0)
  mat = matrix(nrow=n,ncol=u)
  x = x0
  t = t0
  for (i in 1:n) {
    x = stepFun(x,t,times[i]-t,...)
    t = times[i]
    mat[i,] = x
  }
  rownames(mat) = times
  colnames(mat) = names
  mat
}

simSample <- function(n=100, x0, t0=0, deltat, stepFun, ...)
{
  u = length(x0)
  names = names(x0)
  mat = matrix(nrow=n,ncol=u)
  for (i in 1:n) {
    mat[i,] = stepFun(x0,t0,deltat,...)
  }
  colnames(mat) = names
  mat
}


# eof

