gcirc = function(u,ra1,dc1,ra2,dc2) {
  d2r    = pi/180.0e0
  as2r   = pi/648000.0e0
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
  else if(u==2) {
    rarad1 = ra1*d2r
    rarad2 = ra2*d2r
    dcrad1 = dc1*d2r
    dcrad2 = dc2*d2r
  }
  else {
    stop('u must be 0 (radians), 1 ( hours, degrees) or 2 (degrees)')
  }

  deldec2 = (dcrad2-dcrad1)/2.0
  delra2 =  (rarad2-rarad1)/2.0

  sindis = sqrt( sin(deldec2)*sin(deldec2) + 
    cos(dcrad1)*cos(dcrad2)*sin(delra2)*sin(delra2) )
  dis = 2*asin(sindis) 

  if (u!=0)  dis = dis/as2r

  return(dis)
}
