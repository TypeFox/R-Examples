#source('precess.R')

precess_xyz = function(x,y,z,equinox1,equinox2) {

  ra = atan2(y,x)
  del = sqrt(x*x + y*y + z*z)  #magnitude of distance to Sun
  dec = asin(z/del) 
  tmp = precess(ra, dec, equinox1, equinox2, radian=T)
  ra = tmp$ra
  dec = tmp$dec
  equinox1 = tmp$equinox1
  equinox2 = tmp$equinox2
  xunit = cos(ra)*cos(dec)
  yunit = sin(ra)*cos(dec)
  zunit = sin(dec)
  x = xunit * del
  y = yunit * del
  z = zunit * del
  return(list(x=x,y=y,z=z))
}
