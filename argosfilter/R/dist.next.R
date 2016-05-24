.packageName<-"argosfilter"

dist.next<-function(lat, lon)
{
dist.next<-numeric(length(lat))
for (i in 1:(length(lat)-1)) { 
	if (lat[i]==lat[i+1] & lon[i]==lon[i+1]) dist.next[i]<-0 else {
	rlat1=radian(lat[i])
	rlat2=radian(lat[i+1])
	rlon=radian(lon[i+1]-lon[i])
	dist.next[i]<-60*(180/pi)*acos(sin(rlat1)*sin(rlat2)+cos(rlat1)*cos(rlat2)*cos(rlon))
	dist.next[i]<-dist.next[i]*1852 # dist.next is given in m
	}
}
dist.next
}