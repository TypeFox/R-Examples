xyz2ll<-function(x)
{
###  lat here is lattitude (not co-lat)
  RADDEG = 180/pi
  y<-x/sqrt(sum(x*x))
  lat <- asin(y[3])*RADDEG  ;
  lon <- atan2(y[2], y[1])*RADDEG   ;
  return(c(lat,lon))

}
