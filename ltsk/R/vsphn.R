vsphn <-
function(x,sill,range,nugget)
{
  if( range==0) return(rep(sill+nugget,length(x)))
 x1 <- ifelse(x <= range, x, range)
 nugget + sill * ( 3 * x1 / (2*range) - (x1/range)^3/2 )
}
