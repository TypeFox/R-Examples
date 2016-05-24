# function

StepCLE <- function(N,dt=0.01)
{
  S = t(N$Post-N$Pre)
  v = ncol(S)
  sdt = sqrt(dt)
  return(
         function(x0, t0, deltat,...)
         {
           x = x0
           t = t0
           termt = t0+deltat
           repeat {
             h = N$h(x, t, ...)
             dw = rnorm(v,0,sdt)
             dx = S %*% (h*dt + sqrt(h)*dw)
             x = x+as.vector(dx)
             t = t+dt
             if (t > termt)
               return(x)
           }
         }
         )		
}	



# eof

