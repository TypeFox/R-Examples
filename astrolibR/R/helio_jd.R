#source('bprecess.R')
#source('poly.R')
#source('xyz.R')

helio_jd = function(date,ra,dec, B1950=FALSE, time_diff=FALSE) {

  if(!B1950){
    tmp = bprecess(ra,dec)
    ra1 = tmp$ra_1950
    dec1 = tmp$dec_1950
  }
  else {
    ra1 = ra
    dec1 = dec
  }
  
  radeg = 180.0/pi   

  delta_t = (date - 33282.42345905)/36525.0
  epsilon_sec = polyidl( delta_t, c(44.836, -46.8495, -0.00429, 0.00181))
  epsilon = (23.433333e0 + epsilon_sec/3600.0)/radeg
  ra1 = ra1/radeg
  dec1 = dec1/radeg
  tmp=xyz( date)
  x = tmp$x
  y = tmp$y
  z = tmp$z
  
  time = -499.00522*( cos(dec1)*cos(ra1)*x + 
                      (tan(epsilon)*sin(dec1) + cos(dec1)*sin(ra1))*y)
  if(time_diff)
    return(time)
  else 
    return(date + time/86400.0)
}
