vmten <-
function(x,sill,range,nugget)
{
  if( range==0) return(rep(sill+nugget,length(x)))
  xt <- 4.5 * (x / range)
  nugget+ sill *(1-(1+ xt)*exp(-xt))
}
