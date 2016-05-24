# function

StepEulerSPN <- function(N,dt=0.01)
{
  S = t(N$Post-N$Pre)
  return(
         function(x0, t0, deltat,...)
         {
           x = x0
           t = t0
           termt = t0+deltat
           repeat {
             h = N$h(x,t,...)
             x = x+as.vector(S %*% h)*dt
             t = t+dt
             if (t > termt)
               return(x)
           }
         }
         )		
}	



# eof

