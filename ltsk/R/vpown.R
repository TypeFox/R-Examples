vpown <-
function(x,sill,range,nugget)
{
  if( range==0) return(rep(sill+nugget,length(x)))
 nugget+ sill * x^range
}
