#The output is same as that of astrolib mathematically, but IDL prints it transposed.

premat=function( equinox1, equinox2, fk4=F) {

  deg_to_rad = pi/180.0
  sec_to_rad = deg_to_rad/3600.e0
  t = 0.001e0*( equinox2 - equinox1)
  if(!fk4){  
    st = 0.001e0*( equinox1 - 2000.e0)
    a = sec_to_rad * t * (23062.181e0 + st*(139.656e0 +0.0139e0*st) 
      + t*(30.188e0 - 0.344e0*st+17.998e0*t))
    b = sec_to_rad * t * t * (79.280e0 + 0.410e0*st + 0.205e0*t) + a
    c = sec_to_rad * t * (20043.109e0 - st*(85.33e0 + 0.217e0*st) 
      + t*(-42.665e0 - 0.217e0*st -41.833e0*t))
  } else {  
    st = 0.001e0*( equinox1 - 1900.e0)
    a = sec_to_rad * t * (23042.53e0 + st*(139.75e0 +0.06e0*st) 
      + t*(30.23e0 - 0.27e0*st+18.0e0*t))
    b = sec_to_rad * t * t * (79.27e0 + 0.66e0*st + 0.32e0*t) + a
    c = sec_to_rad * t * (20046.85e0 - st*(85.33e0 + 0.37e0*st) 
      + t*(-42.67e0 - 0.37e0*st -41.8e0*t))
  } 
  sina = sin(a)
  sinb = sin(b)
  sinc = sin(c)
  cosa = cos(a)
  cosb = cos(b)
  cosc = cos(c)
  r = matrix(
  c( cosa*cosb*cosc-sina*sinb, sina*cosb+cosa*sinb*cosc,  cosa*sinc,
  -cosa*sinb-sina*cosb*cosc, cosa*cosb-sina*sinb*cosc, -sina*sinc,
  -cosb*sinc, -sinb*sinc, cosc),3,3)
  return(r)
}
