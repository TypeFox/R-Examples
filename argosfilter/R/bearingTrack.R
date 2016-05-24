.packageName<-"argosfilter"

`bearingTrack` <-
function(lat, lon)
{
bearingTrack<-numeric(length(lat)-1)
for (k in 1:(length(lat)-1)){
	lat1<-lat[k]
	lat2<-lat[k+1]
	lon1<-lon[k]
	lon2<-lon[k+1]
	bearingTrack[k]<-bearing(lat1,lat2,lon1,lon2)
	}
bearingTrack
}

