#This is a wrapper for the georeference function for multiple GPS locations within a Dead Reckoning Track
GeoRef<-function(drdata, gpsfdata)
{
	timer=proc.time()[3]
	
	
	matched<-match(gpsfdata$DateTime,drdata$DateTime)
	matched<-matched[!is.na(matched)]
	matched2<-match(drdata$DateTime,gpsfdata$DateTime)
	matched2<-unique(matched2[!is.na(matched2)])
	Georeferenced<-as.data.frame(matrix(-666,max(matched),10))
	colnames(Georeferenced)<-c("Distance","LatRad","LongRad","Latitude","Longitude","Depth","DRCalcSpd","NewX","NewY","Bearing")
	DateTime<-character(length=max(matched))
	print(c("1",proc.time()[3]-timer))
	for(i in 1:(length(matched)-1))
	{
		startx<-drdata[matched[i]:matched[i+1],]$Xdim[1]
		drdata[matched[i]:(matched[i+1]-1),]$Xdim<-drdata[matched[i]:(matched[i+1]-1),]$Xdim-startx
		starty<-drdata[matched[i]:matched[i+1],]$Ydim[1]
		drdata[matched[i]:(matched[i+1]-1),]$Ydim<-drdata[matched[i]:(matched[i+1]-1),]$Ydim-starty
		LL<-dim(drdata[matched[i]:(matched[i+1]-1),])[1]
		tempGR<-GeoReference(drdata[matched[i]:(matched[i+1]-1),],gpsfdata[matched2[i]:matched2[i+1],])[1:LL,]
		Georeferenced[matched[i]:(matched[i+1]-1),]<-tempGR[,2:11]
		DateTime[matched[i]:(matched[i+1]-1)]<-tempGR[,1]
		print(c(i,proc.time()[3]-timer))
	}
	Georeferenced<-cbind(DateTime,Georeferenced,stringsAsFactors=F)
	Georeferenced<-Georeferenced[Georeferenced$LatRad!=(-666),]
	return(Georeferenced=Georeferenced)
}
