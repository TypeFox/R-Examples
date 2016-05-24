vmte <-
function(x,sill,range)
{
  if( range==0) return(rep(sill,length(x)))
  xt <- 4.5 * (x / range)
  sill *(1-(1+ xt)*exp(-xt))
}
