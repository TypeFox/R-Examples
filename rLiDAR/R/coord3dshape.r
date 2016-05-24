# Author: Remko Duursma. Uses code by Remko Duursma in Maeswrap package, see Plotstand.
coord3dshape <- function(crownshape=c("cone","ellipsoid","halfellipsoid","paraboloid","cylinder"),
                         nz=5, nalpha=5, CL=1, CW=1, x0=0, y0=0, z0=0
){
  
  crownshape <- match.arg(crownshape)
  
  z <- rep(seq(0,1,length=nz),each=nalpha)
  angs <- rep(seq(0,2*pi, length=nalpha),nz)
  
  if(crownshape == "cone")distfun <- (1-z)
  if(crownshape == "ellipsoid")distfun <- sqrt(1 - ((z-1/2)^2)/((1/2)^2))
  if(crownshape == "halfellipsoid")distfun <- sqrt(1 - z**2)
  if(crownshape == "paraboloid")distfun <- sqrt(1-z)
  if(crownshape == "cylinder")distfun <- 1
  
  
  r <- CW/2
  x <- x0 + r*distfun*cos(angs)
  y <- y0 + r*distfun*sin(angs)
  z <- z0 + z*CL
  
  keep <- !duplicated(cbind(x,y,z))
  x <- x[keep]
  y <- y[keep]
  z <- z[keep]
  return(matrix(cbind(x,y,z),ncol=3))
}



