ABEL1Q <-
function(T,C)
{
   ALPHA <- 4 * (C[1]^2)/((1-C[1])^4) 
   ALPHA <- 1.3221*((ALPHA*T)^.2)
return(ALPHA)
}
