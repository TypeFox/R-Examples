.packageName<-"argosfilter"

`bearing` <-
function(lat1, lat2, lon1, lon2)
{
if (lon1==lon2 & lat1==lat2) bearing<-NA else{
	if (lon1==lon2) {if (lat1<=lat2) bearing<-0 else bearing<-180}
	if (lon1!=lon2) {
		d<-distance_m(lat1,lon1,lat2,lon2)/1852
		rlat1=radian(lat1)
		rlat2=radian(lat2)
		rlon=radian(lon2-lon1)
		bearing_rad<-acos((sin(rlat2)-sin(rlat1)*cos(radian(d/60)))/(sin(radian(d/60))*cos(rlat1)))
		bearing=bearing_rad*180/pi
		if(sin(rlon)<0) bearing=360-bearing
		if(bearing>180) bearing=bearing-360
		}	
	}
bearing
}

