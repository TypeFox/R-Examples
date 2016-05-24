legend.bubbles <-
function (x, y=NULL, z=NULL, nleg=NULL, digits=NULL, pch, 
z.cex.range=NULL, x.intersp=1, y.intersp=1, bg="white", ...)

{
  if ( is.null ( x )) stop ("x is missing")
  if ( is.null ( z )) stop ("z is missing")
  if ( is.null ( z.cex.range )) stop ("z.cex.range is missing")
  if ( is.null ( nleg )) { nleg<-3 }
  if ( is.null ( digits )) { digits<-1 }

  n <- (( 1 ) / ( nleg-1 ))
  probs <- seq ( 0 , 1 , n )
  z <- as.numeric ( round ( quantile ( z ,  probs = probs ) , digits = digits ))

  cex.bubbles <- as.numeric ( TT.str ( z , z.cex.range [ 1 ] , z.cex.range [ 2 ] ))
  
  legend ( x = x , y = y , legend = z , pch = pch , pt.cex = cex.bubbles , x.intersp = x.intersp , y.intersp = y.intersp , bg=bg,...)
}
