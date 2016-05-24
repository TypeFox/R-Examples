vexp <-
function(x,sill,range)
{
  if( range==0) return(rep(sill,length(x)))
  sill *( 1 - exp(-3*x/range))
}
