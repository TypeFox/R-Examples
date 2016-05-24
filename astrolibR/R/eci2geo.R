#source('ct2lst.R')

eci2geo = function(eci_xyz,jdtim) {
  if(!is.matrix(eci_xyz)) eci_xyz = matrix(eci_xyz,1)
  re=6378.137     # earth's equatorial radius, in km
  coord=eci_xyz
  jdtime= jdtim
  theta=atan2(coord[,2],coord[,1]) # azimuth       
  gst = ct2lst(0,0,jdtime)

  angle_sid=gst*2*pi/24 # sidereal angle
  lon= (theta - angle_sid ) %% (2*pi) #longitude      
  r=sqrt(coord[,1]^2+coord[,2]^2)
  lat=atan2(coord[,3],r) # latitude
  alt=r/cos(lat) - re # altitude 
  lat=lat*180./pi # to convert from radians into degrees...
  lon=lon*180./pi
  ss=(lon<0) 
  lon[ss]=lon[ss]+360
  
  return (c(lat,lon,alt))
}
  
