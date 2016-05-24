### cart2clock.R

cart2clock <- function(x, y = NULL, circle)
  ## Author: Rene Locher
  ## Version: 2005-10-05
  ##
  ## converts clock coordinates to cartesian coordinates
  ## Attention: phi corresponds to the angle
  ## between y and the direction wanted, measured clockwise.
  ## This corresponds to an angle pi/2-phi in ordinary polar coordinates
{
  xy <- getXY(x, y, unidim.allowed=FALSE)
  x <- xy$x
  y <- xy$y
  phi <-(atan(x/y)/2/pi*circle + ifelse(y>=0,circle,1.5*circle))%%circle
  return(data.frame(rho=sqrt(x*x+y*y),phi=phi))
} ## cart2clock
