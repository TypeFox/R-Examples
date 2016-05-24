normalize <- function( x, base=2 ) {
##  	 ##*   return  a and e so that  x  = a*base^e,   abs( a ) in [1, base)
  e <- trunc(log(abs(x),base=base))
  if (e == Inf) e <- 0
  a <- x/base^e
  if ( abs(a) < 1.0 ) { e <- e-1.0;   a <- a*base
  } else { if (abs(a) > base ){ e <- e+1.0; a <- a/base }
  }
  return(c(a,e))
} ## end  normalize
