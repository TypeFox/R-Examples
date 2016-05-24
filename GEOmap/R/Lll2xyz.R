Lll2xyz<-function(lat, lon)
{
 ### lat here is lattitude (not co-lat)
  DEGRAD = pi/180
  z = sin(DEGRAD*lat);
  y = cos(DEGRAD*lat)*sin(DEGRAD*lon);
  x = cos(DEGRAD*lat)*cos(DEGRAD*lon);

  return(list(x=x,y=y,z=z))
}
