#Function to do the 2d interpolation (on the log scale)

'interplcd' <-
  function(lcd, gridlen=100 )
{
  x <- lcd$x[,1]
  y <- lcd$x[,2]
  z <- lcd$logMLE
  
  xo <- seq(min(x), max(x), length=gridlen)
  yo <- seq(min(y), max(y), length=gridlen)

  chull <- lcd$triang
  nsimplices <- nrow(chull)

  grid <- as.matrix( expand.grid( xo, yo ) )
  isout <- apply( lcd$outnorm %*% ( t( grid ) - lcd$midpoint ) -
                 lcd$outdist, 2, max ) > 0
  g <- apply( lcd$bunique %*% t( grid )  - lcd$betaunique, 2, min ) +
  ifelse( isout, -Inf, 0 )
  g <- matrix( g, nrow=gridlen, ncol=gridlen )
  
  return(list(x=xo,y=yo,z=g))
}
