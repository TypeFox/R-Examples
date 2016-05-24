.packageName<-"argosfilter"

`distance` <-
function(lat1, lat2, lon1, lon2)
{
if (lat1==lat2 & lon1==lon2) distance<-0 else {
	rlat1=radian(lat1)
	rlat2=radian(lat2)
	rlon=radian(lon2-lon1)
	distance<-60*(180/pi)*acos(sin(rlat1)*sin(rlat2)+cos(rlat1)*cos(rlat2)*cos(rlon))
	distance<-distance*1852/1000 #distance is given in km
	}
distance
}

