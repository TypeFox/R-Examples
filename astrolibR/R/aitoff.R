aitoff = function(l,b) {
  radeg = 180/pi
  sa = l
  x180 = which (sa>180.0)
  sa[x180]  = sa[x180] - 360.
  alpha2 = sa/(2*radeg)
  delta = b/radeg   
  r2 = sqrt(2)    
  f = 2*r2/pi     
  cdec = cos(delta)    
  denom =sqrt(1. + cdec*cos(alpha2))
  x = cdec*sin(alpha2)*2.*r2/denom
  y = sin(delta)*r2/denom
  x = x*radeg/f
  y = y*radeg/f
  return (list(x=x,y=y))
}
