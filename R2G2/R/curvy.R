curvy <-
function(f, startDD, stopDD){  
  radian = 180 / pi

  lat1 = startDD[2] / radian
  lon1 = startDD[1] / radian 
  lat2 = stopDD[2] / radian
  lon2 = stopDD[1] / radian

  d = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2))
  A = sin((1 - f) * d) / sin(d)
  B = sin(f * d) / sin(d)
  x = as.numeric(A * cos(lat1) * cos(lon1) +  B * cos(lat2) * cos(lon2))
  y = as.numeric(A * cos(lat1) * sin(lon1) +  B * cos(lat2) * sin(lon2))
  z = as.numeric(A * sin(lat1) +  B * sin(lat2))
  lat=atan2(z, sqrt(x^2 + y^2))
  lon=atan2(y, x)
  rev(c(lat * 180/pi, lon * 180/pi))
  }

