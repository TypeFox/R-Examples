Lxyz2ll<-function(X)
{
###  lat here is lattitude (not co-lat)
  x = cbind(X$x, X$y, X$z)
  r = sqrt( X$x^2+ X$y^2 + X$z^2 )
  y = cbind(x[,1]/r, x[,2]/r, x[,3]/r)
  
  lat <- asin(y[,3])*180/pi  ;
  lon <- atan2(y[,2], y[,1])*180/pi   ;

  return(list(lat=lat,lon=lon))
}
