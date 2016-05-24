
require(deSolve)

StepODE <- function(RHSfun)
{
  return(
         function(x0, t0, deltat, ...)
         {
           termt = t0+deltat
           res = deSolve::ode(x0,times=c(t0,termt),func=function(t, x, ...) list(RHSfun(x, t, ...)),...)
           res = res[2,]
           res[2:length(res)]
         }
         )		
}	



# eof

