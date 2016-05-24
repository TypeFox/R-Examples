
StepPTS <- function(N,dt=0.01)
{
  S = t(N$Post-N$Pre)
  v = ncol(S)
  return(
         function(x0,t0,deltat,...)
         {
           x = x0
           t = t0
           termt = t0+deltat
           repeat {
             h = N$h(x, t, ...)
             r = rpois(v,h*dt)
             x = x+as.vector(S %*% r)
             t = t+dt
             if (t > termt)
               return(x)
           }
         }
         )		
}	



# eof

