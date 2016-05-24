"Ucutoff" <-
function(p, q1, q2, VarU)
{
   ncp = (q2-q1)^2 / 4 / VarU
   return( sqrt(qchisq(p, 1, ncp)) )
}

