#Convert between lat, long in degree and northing/easting in km.

DegToKM <- function(gpsformat)
{
	#Function to further convert the formatted GPS points into the coordinates
	#Unites in meters.
	K <- nrow(gpsformat)
	Xdiff<- gpsformat$DistanceKm[-1]*sin(gpsformat$BearingRad[-K])
	Ydiff<- gpsformat$DistanceKm[-1]*cos(gpsformat$BearingRad[-K])
	
	gpsX <- c(0, cumsum(Xdiff))
	gpsY <- c(0, cumsum(Ydiff))
	return(data.frame(DateTime=gpsformat$DateTime, Easting=gpsX, Northing=gpsY))
}

KMToDeg <- function(cPath, iniDeg)
{
	#Use Brian's function to convert the points back to Latitude and Longitude.
	#gdr, first column X (easting), next column Y(northing). Unit in KM!!!!!
	#iniDeg, longitude first and latitude second.
	gdr <- cPath*1000
	T <- nrow(gdr)
	gdr <- as.matrix(gdr)
	ini <- as.numeric(iniDeg)
	diffMx <- gdr[-1, ] - gdr[-T, ]
	Distance <- sqrt(diffMx[, 1]^2 + diffMx[, 2]^2)
    Bering <- atan2(diffMx[,1], diffMx[,2])
	degMx <- matrix(NA, T, 2)
	initLat <- ini[2]
	initLong <- ini[1]
	degMx[1, ] <- c(initLat, initLong)
	for (i in 2:T)
	{
		tLatRad <- CalcLatitude(initLat, Distance[i-1], Bering[i-1])
		tLonRad <-CalcLongitude(initLat, tLatRad, initLong, 
			Distance[i-1], Bering[i-1])
		degMx[i, ] <- c(tLatRad, tLonRad)*180/pi
		initLat <- degMx[i, 1]
		initLong <- degMx[i, 2]
	}
	return(data.frame(Longitude=degMx[,2], Latitude=degMx[,1]))
}