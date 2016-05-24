GeoReference <-function(drdata, gpsdata)
{
	###########drdata first entry must line up with the time stamp of the first engty for GPS data#####################
	#Based on Wilson et al 2007  All at sea with Animal tracks...
	#Requires "Navigation Functions"
	#timer=proc.time()[3] # used to determine slow parts of the program and show some output while it is running to indicate where a failure occurs.
	drnrow <- nrow(drdata)
	
	GPSBearing=gpsdata[1,6]
	GPSDistance <- gpsdata[2,8]
	GPSXLong <- sin(GPSBearing)*GPSDistance*1000
	GPSYLat <- cos(GPSBearing)*GPSDistance*1000
	initLat <- gpsdata$Latitude[1] #gpsdata Lat Degrees
	initLong <- gpsdata$Longitude[1] #gpsdata Long Degrees
	#print(c("2",proc.time()[3]-timer))
	FinalDRdistance <- sqrt(drdata$Xdim[drnrow]^2+drdata$Ydim[drnrow]^2)
	FinalDRBearing<-atan2(drdata$Xdim[drnrow],drdata$Ydim[drnrow])
	FinalDRLat<-CalcLatitude(initLat,FinalDRdistance,FinalDRBearing)
	FinalDRLong<-CalcLongitude(initLat,FinalDRLat,initLong,FinalDRdistance,FinalDRBearing)
	#print(c("3",proc.time()[3]-timer))
	AdjXLong=(GPSXLong-drdata$Xdim[drnrow])/(sum(drdata$Speed)) #drnrow-1
	AdjYLat=(GPSYLat-drdata$Ydim[drnrow])/(sum(drdata$Speed)) #drnrow-1
	#print(c("4",proc.time()[3]-timer))
	#rownum=c(0,seq(drnrow-1))
	Spdsum=cumsum(drdata$Speed)
	#print(c("5.1",proc.time()[3]-timer))
	#GeoRefLatLong <- matrix(0,drnrow,11,dimnames=list(1:drnrow,c("DateTime","Distance","LatRad","LongRad","Latitude","Longitude","Depth","DRCalcSpd","NewX","NewY","Bearing")))
	#print(c("5.2",proc.time()[3]-timer))
	NewX<-drdata$Xdim+AdjXLong*Spdsum #rownum
	NewY<-drdata$Ydim+AdjYLat*Spdsum #rownum
	#print(c("5.3",proc.time()[3]-timer))
	Distance<-sqrt(NewX[2:drnrow]^2+NewY[2:drnrow]^2)
	Bearing<-atan2(NewX[2:drnrow],NewY[2:drnrow])
	#GeoRefLatLong[,9] <- NewX#format(NewX, digits=10, drop0trailing=FALSE )#NewX
	#rm(NewX)
	#GeoRefLatLong[,10] <- NewY#format(NewY, digits=10, drop0trailing=FALSE )#NewY
	#rm(NewY)
	#print(c("5.4",proc.time()[3]-timer))
	LatRad<-CalcLatitude(initLat,Distance,Bearing)
	LongRad<-CalcLongitude(initLat,LatRad,initLong,Distance,Bearing)
	#print(c("5.5",proc.time()[3]-timer))
	LatRad<-c(initLat*pi/180,LatRad)
	LongRad<-c(initLong*pi/180,LongRad)
	#print(c("5.6",proc.time()[3]-timer))
	Bearing<-c(Bearing,0)
	#GeoRefLatLong[,11] <- Bearing#format(Bearing, digits=10, drop0trailing=FALSE )#Bearing
	#rm(Bearing)
	Distance<-c(0,Distance)
	#print(c("6",proc.time()[3]-timer))
	
	Latitude<-LatRad/(pi)*180
	Longitude<-LongRad/(pi)*180
	#print(c("7",proc.time()[3]-timer))
	
	#XYZ_Distance=c(0,sqrt(drdata$Speed[2:drnrow]^2+(Distance[2:drnrow]-Distance[1:drnrow-1])^2))
	#GeoRefLatLong[,9] <- XYZ_Distance#format(XYZ_Distance, digits=10, drop0trailing=FALSE )#XYZ_Distance
	#rm(XYZ_Distance)
	#GeoRefLatLong[,2] <- Distance#format(Distance, digits=10, drop0trailing=FALSE )#Distance
	#rm(Distance)
	#GeoRefLatLong[,3] <- LatRad#format(LatRad, digits=10, drop0trailing=FALSE)
	#rm(LatRad)
	#GeoRefLatLong[,4] <- LongRad#format(Longitude, digits=10, drop0trailing=FALSE)
	#rm(LongRad)
	#GeoRefLatLong[,5] <- Latitude#format(Latitude, digits=10, drop0trailing=FALSE)
	#rm(Latitude)
	#GeoRefLatLong[,6] <- Longitude#format(LongDeg, digits=10, drop0trailing=FALSE)
	#rm(Longitude)
	#GeoRefLatLong[,7] <- drdata$Depth#format(drdata$Depth, digits=5, drop0trailing=FALSE )#Depth
	#GeoRefLatLong[,8] <- format(drdata[,5], digits=10, drop0trailing=FALSE )#ODBA
	#GeoRefLatLong[,8] <- drdata$Speed#format(drdata$Speed, digits=10, drop0trailing=FALSE, trim=FALSE)#DRCalcSpeed
	#print(c("8",proc.time()[3]-timer))
	#Format the TimeDate stamp	
	#GeoRefLatLong <- as.data.frame(GeoRefLatLong)
	#print(c("9",proc.time()[3]-timer))
	#GeoRefLatLong[,1] <- as.character(drdata$DateTime)
	GeoRefLatLong<- data.frame(as.character(drdata$DateTime),Distance,LatRad,LongRad,Latitude,Longitude,drdata$Depth,drdata$Speed,NewX,NewY,Bearing)
	names(GeoRefLatLong)<-c("DateTime","Distance","LatRad","LongRad","Latitude","Longitude","Depth","DRCalcSpd","NewX","NewY","Bearing")
	#print(c("10",proc.time()[3]-timer))
	#print(c("1",proc.time()[3]-timer)) # Testing and failure information
	return(GeoRefLatLong=GeoRefLatLong)#list(AdjXLong=AdjXLong,AdjYLat=AdjYLat))#for Squish factors
}
