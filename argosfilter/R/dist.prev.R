.packageName<-"argosfilter"

dist.prev<-function(lat, lon)
{
dist.prev<-numeric(length(lat))
for (i in 2:length(lat)) { 
	if (lat[i]==lat[i-1] & lon[i]==lon[i-1]) dist.prev[i]<-0 else {
	rlat1=radian(lat[i])
	rlat2=radian(lat[i-1])
	rlon=radian(lon[i-1]-lon[i])
	dist.prev[i]<-60*(180/pi)*acos(sin(rlat1)*sin(rlat2)+cos(rlat1)*cos(rlat2)*cos(rlon))
	dist.prev[i]<-dist.prev[i]*1852 # dist.prev is given in m
	}
}
dist.prev
}

