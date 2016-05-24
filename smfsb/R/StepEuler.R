# function

StepEuler <- function(RHSfun,dt=0.01)
{
  return(
         function(x0, t0, deltat,...)
         {
           x = x0
           t = t0
           termt = t0+deltat
           repeat {
             x = x+RHSfun(x, t, ...)*dt
             t = t+dt
             if (t > termt)
               return(x)
           }
         }
         )		
}	



# eof

