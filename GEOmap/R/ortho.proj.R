ortho.proj<-function(lat, lon, lon0, lat1, R)
  {

    if(missing(R))  R = 6378.2064
    
    
    lam = (lon-lon0)*pi/180
    phi = lat*pi/180
    phi1 =  lat1*pi/180
    x = R*cos(phi)*sin(lam)
    y = R*(cos(phi1)*sin(phi)  - sin(phi1)*cos(phi)*cos(lam))

    cosc = sin(phi1)*sin(phi) +cos(phi1)*cos(phi)*cos(lam)

    return(list(x=x, y=y, cosc=cosc))

  }
