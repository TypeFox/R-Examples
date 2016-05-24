longlat.to.km <-
function(long1,lat1,long2,lat2){
# compute surface distances in k between two points with geographical coordinates 
# long1,lat1 and long2,lat2 
# not used in the first version of package etasFLP
radius=6371.3
long1=long1*pi/180
long2=long2*pi/180
lat1=lat1*pi/180
lat2=lat2*pi/180

return(radius*acos(cos(long1-long2)*cos(lat1)*cos(lat2)+sin(lat1)*sin(lat2)))
}
