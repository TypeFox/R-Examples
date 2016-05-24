vpow <-
function(x,sill,range)
{
  if( range==0) return(rep(sill,length(x)))
  sill * x^range
}
