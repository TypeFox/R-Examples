"SE.Meas" <-
function (s, rxx)
{
sem <- s*sqrt(1-rxx)
return(sem)
 }

