`new.lat.long` <-
function(long,lat,bearing,distance) {
  rad <- pi/180
  a1<-lat*rad
  a2<-long*rad
  tc<-bearing*rad
  d<-distance/6378.1
  nlat<-asin(sin(a1)*cos(d)+cos(a1)*sin(d)*cos(tc))
  dlon<-atan2(sin(tc)*sin(d)*cos(a1),cos(d)-sin(a1)*sin(nlat))
  nlon<-((a2+dlon+pi)%%(2*pi))-pi
  npts<-c(nlat/rad,nlon/rad)
  return(npts)
}

