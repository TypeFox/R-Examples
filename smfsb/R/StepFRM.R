
StepFRM <- function(N)
{
  S = t(N$Post-N$Pre)
  v = ncol(S)
  return(
         function(x0, t0, deltat, ...)
         {
           x = x0
           t = t0
           termt = t0+deltat
           repeat {
             h = N$h(x, t, ...)
             pu = rexp(v,h)
             j = which.min(pu)
             t = t+pu[j]
             if (t >= termt)
               return(x)
             x = x+S[,j]
           }
         }
         )		
}	


# eof

