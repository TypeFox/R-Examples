vgaun <-
function(x,sill,range,nugget)
{
  if( range==0) return(rep(sill+nugget,length(x)))
  nugget + sill* (1 - exp(-3*(x/range)^2))
}
