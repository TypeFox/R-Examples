precess = function( ra, dec, equinox1, equinox2,
                    fk4=F,radian=F) {

  deg_to_rad = pi/180

  if(!radian ){
    ra_rad = ra*deg_to_rad     #convert to double
                                        #precision if not already
    dec_rad = dec*deg_to_rad 
  }
  else {
    ra_rad= ra
    dec_rad = dec
  }
  a = cos( dec_rad )  
    x = rbind(a*cos(ra_rad),
              a*sin(ra_rad),
              sin(dec_rad))

  sec_to_rad = deg_to_rad/3600.e0
  r = premat(equinox1, equinox2, fk4 = fk4)
  cat('precess:','equinox1=',equinox1,'\n')
  cat('precess:','equinox2=',equinox2,'\n')

  x2 = r%*%x      #rotate to get output direction cosines

  cat('precess:','x2=',x2,'\n')

    ra_rad = atan2(x2[2,],x2[1,])
    dec_rad = asin(x2[3,])

  if(!radian ){
    ra = ra_rad/deg_to_rad
    ra = ra + (ra<0.)*360            #ra between 0 and 360 degrees
    dec = dec_rad/deg_to_rad
  } 
  else {
    ra = ra_rad
    dec = dec_rad
    ra = ra + (ra<0.)*2*pi
  }
  return(list(ra =ra, dec=dec))
}
