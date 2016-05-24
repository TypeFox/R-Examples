`wiggle.env` <-
function(x,  y )
{
  ###  givn an x,y signal, get the smoothed peaks
  why = y-mean(y)	
  aa = peaks(why, span=3, do.pad = TRUE)
  w = aa & why>0

  smsp = smooth.spline( cbind(x[w], why[w])  , spar = 0.2)

  return( smsp )
}

