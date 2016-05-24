GPStable <-
function(rawdata)
{
	#Format of raw data should be "DateTime" "Latitude" "Longitude"
	Nrow <- nrow(rawdata)
	
	gpsdata <- matrix(0,Nrow,8)
		colnames(gpsdata)<-c("DateTime","Latitude","Longitude","LatRad","LongRad","BearingRad","BearingDeg","DistanceKm")
	gpsdata <-data.frame(gpsdata)
	
	#gpsdata[,1] <- 0:(Nrow-1)
	if('DateTime' %in% colnames(rawdata))
	{
		gpsdata[,1]<-rawdata$DateTime
	}else{
		gpsdata[,1]<-paste(as.character(rawdata$Date),as.character(rawdata$Time)) #Combine date and time
	}
	gpsdata[,2] <- rawdata$Latitude
	gpsdata[,3] <- rawdata$Longitude
	#Turn Degrees into Radians
	gpsdata[,4]<- gpsdata[,2]/360*2*pi
	gpsdata[,5]<-gpsdata[,3]/360*2*pi
	#This is for the bering and requires the final one to be 0
	for(i in 1:(Nrow-1))
	{
		initialLong<-gpsdata[i,3]
		finalLong<-gpsdata[i+1,3]
		initialLat<-gpsdata[i,2]
		finalLat<-gpsdata[i+1,2]
		#Radians
		gpsdata[i,6]<-CalcBearing(initialLat,initialLong,finalLat,finalLong)
		#Degrees
		gpsdata[i,7]<-ifelse(gpsdata[i,6]/(2*pi)*360>0,gpsdata[i,6]/(2*pi)*360,360+gpsdata[i,6]/(2*pi)*360)
	}
	#Distance using the Spherical Law of Cosines
	for(i in 2:Nrow)
	{
		initialLong<-gpsdata[i-1,3]
		finalLong<-gpsdata[i,3]
		initialLat<-gpsdata[i-1,2]
		finalLat<-gpsdata[i,2]
		gpsdata[i,8]<-CalcDistance(initialLat, initialLong, finalLat, finalLong)
	}
	
	return(gpsdata=gpsdata)
}
