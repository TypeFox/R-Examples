posang = function(u,ra1,dc1,ra2,dc2) {

  if (all(ra1==ra2) && all(dc1==dc2))
    stop(paste("positions are equal:", ra1, dc1))


  d2r    = pi/180.0e0
  h2r    = pi/12.0e0

  if(u==0) {
    rarad1 = ra1
    rarad2 = ra2
    dcrad1 = dc1
    dcrad2 = dc2
  }
  else if(u==1) {
    rarad1 = ra1*h2r
    rarad2 = ra2*h2r
    dcrad1 = dc1*d2r
    dcrad2 = dc2*d2r
  }
  else {
    stop('u must be 0 for radians or 1 for hours, degrees, arcsec')
  }
  radif  = rarad2-rarad1
  angle  = atan2(sin(radif),
                cos(dcrad1)*tan(dcrad2)-sin(dcrad1)*cos(radif))
  if (u!=0) angle = angle/d2r  

  return(angle)
}             
