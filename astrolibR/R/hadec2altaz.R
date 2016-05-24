hadec2altaz = function( ha, dec, lat, ws=F) {

  d2r = pi/180.
  sh = sin(ha*d2r) 
  ch = cos(ha*d2r)
  sd = sin(dec*d2r) 
  cd = cos(dec*d2r)
  sl = sin(lat*d2r) 
  cl = cos(lat*d2r)
  x = - ch * cd * sl + sd * cl
  y = - sh * cd
  z = ch * cd * cl + sd * sl
  r = sqrt(x^2 + y^2)
  az = atan2(y,x) /d2r
  alt = atan2(z,r) / d2r
  w = which(az<0)
  az[w] = az[w] + 360.
  if(ws) az = (az + 180.) %% 360.
  return(list(alt=alt,az=az))
}
