bandwidth.parameter <-
function(p,n)
{
 return(round((4/((p+2)*n))^(1/(p+4)),3))
}
