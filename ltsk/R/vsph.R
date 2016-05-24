vsph <-
function(x,sill,range)
{
  if( range==0) return(rep(sill,length(x)))
 x1 <- ifelse(x <= range, x, range)
 sill * ( 3 * x1 / (2*range) - (x1/range)^3/2 )
}
