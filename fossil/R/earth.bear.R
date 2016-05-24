`earth.bear` <-
function(long1, lat1, long2, lat2) {
  rad <- pi/180
  a1 <- lat1*rad
  a2 <- long1*rad
  b1 <- lat2*rad
  b2 <- long2*rad
  dlon <- b2-a2
  bear <- atan2(sin(dlon)*cos(b1),cos(a1)*sin(b1)-sin(a1)*cos(b1)*cos(dlon))
  deg <- (bear%%(2*pi))*(180/pi)
  return(deg)
}

