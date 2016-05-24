"Wcutoff" <-
function(p, eps1, eps2, VarW)
{
   ncp <- (eps1+eps2)^2/4/VarW
   return(sqrt(qchisq(p, 1, ncp)))
}

