along.great<-function( phi1, lam0,  c ,  Az  )
  {
#########         given a point (lat lon) radians
########         a distance in radians (c) and
########         an azimuthal direction Az
########         calculate the point (phi, lam) which lies angle radians away from (lat, lon)
    

    phi =  asin( sin(phi1)*cos(c) + cos(phi1)*sin(c)*cos(Az));
    tem = atan2( sin(c)*sin(Az), (cos(phi1)*cos(c) - sin(phi1)*sin(c)*cos(Az)));
    lam = lam0 + tem;

    return(list(phi=phi, lam=lam))
  }
