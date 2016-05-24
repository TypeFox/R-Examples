getgreatarc<-function( lat1,  lon1,  lat2,  lon2,  num=10)
  {
    if(missing(num)) { num=10 }
    ##  get an array of points along the great circle between two points
    RADDEG = 180/pi
    DEGRAD = pi/180
    daz = distaz(lat1,  lon1,  lat2,  lon2)
    dis = daz$del
    Az = daz$az
    dang = dis/num;
    angle = seq(from=0, to=dis*pi/180, length=num)  
    AA = along.great( lat1*DEGRAD,  lon1*DEGRAD, angle, Az*DEGRAD);
    phi2 = RADDEG * AA$phi;
    lam2 = RADDEG * AA$lam;
    return(list( lat=phi2, lon=lam2, del=dis ) )

    
  }


