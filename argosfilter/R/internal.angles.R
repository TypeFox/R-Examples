.packageName<-"argosfilter"

internal.angles<-function(lat, lon)
{

# calculate bearing to the next position
bn<-bearingTrack(lat,lon)
bn[1:10]
bn<-c(bn,0)
# calculate bearing to the previous position
bp2<-bearingTrack(rev(lat),rev(lon))
bp<-rev(bp2)
bp<-c(0,bp)
bp[1:10]
# calculate the angle to the previous and next position -----
ang<-numeric(length(lat))
for (i in 2:(length(lat)-1)) {
		if (lat[i]==lat[i+1]&lon[i]==lon[i+1] | lat[i]==lat[i-1]&lon[i]==lon[i-1]) ang[i]<-180 else {
		ang[i]<-abs(bp[i]-bn[i])
			if (ang[i]>180) ang[i]<-(360-ang[i]) else ang[i]<-ang[i] 
			}
		}

ang[1:10]
ang
}
