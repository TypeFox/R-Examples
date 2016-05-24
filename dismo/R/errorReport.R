# Download geographic data and return as R object
# Author: Jacob van Etten
# License GPL3
# Version 0.1
# October 2008


.errorReport <- function(xyxy,threshold) 
{	
	distn <- matrix(ncol=16,nrow=length(xyxy[,1]))
	for(i in 1: length(xyxy[,1])) {
		distn[i,] <- .errorDist(xyxy[i,], threshold)
	}
	errorNames <- c("Imprecision","Lonlat swap","Sign latitude","Sign longitude","Sign latitude and longitude","Lonlat swap and sign latitude","Lonlat swap and sign longitude","Lonlat swap and both signs","Wrong longitude","Wrong latitude","Lonlat swap and wrong longitude","Lonlat swap and wrong latitude","Sign latitude and wrong longitude","Sign longitude and wrong latitude","Sign latitude, lonlat swap and wrong longitude","Sign longitude, lonlat swap and wrong latitude")
	classify <- function(x) {
		if (any(is.na(x))){out <- NA}
		else {out <- errorNames[min(which(x==min(x)))]}
		return(out)
	}
	
	errorCategory <- apply(distn,1,classify)
	
	classify2 <- function(x) {
		if (any(is.na(x))){out <- NA}
		else {out <- min(which(x==min(x)))}
		return(out)
	}
	index <- apply(distn,1,classify2)
	index <- cbind(c(1:length(xyxy[,1])),index)
	
	dist2 <- matrix(ncol=16,nrow=length(xyxy[,1]))
	for(i in 1: length(xyxy[,1])) 	{
		dist2[i,] <- .errorDist(xyxy[i,],0)
	}
	
	errorDistance <- dist2[index]/1000
	result <- cbind(errorCategory,errorDistance)
	return(result)
}


.errorDist <- function(pp,a){
		x1 <- pp[1]
		y1 <- pp[2]
		x2 <- pp[3]
		y2 <- pp[4]

		if(any(is.na(c(x1, y1, x2, y2)))){
			out <- NA
		} else {	
			difference <- pointDistance(c(x1, y1), c(x2, y2), type='GreatCircle') 
			xy.exchange <- pointDistance(c(y1, x1), c(x2, y2), type='GreatCircle') 
			signch.lat <- pointDistance(c(x1, -y1), c(x2, y2), type='GreatCircle') 
			signch.lon <- pointDistance(c(-x1, y1), c(x2, y2), type='GreatCircle')
			signch.latlon <- pointDistance(c(-y1, -x1), c(x2, y2), type='GreatCircle')
			
			xy.signch.lat <- pointDistance(c(y1, -x1), c(x2, y2), type='GreatCircle')
			xy.signch.lon <- pointDistance(c(-y1, x1), c(x2, y2), type='GreatCircle')
			xy.signch.latlon <- pointDistance(c(-y1, -x1), c(x2, y2), type='GreatCircle')
			
			wrong.lon <- pointDistance(c(0, y1), c(a, y2), type='GreatCircle')
			wrong.lat <- pointDistance(c(x1, 0), c(x2, a), type='GreatCircle')
			
			xy.wrong.lon <- pointDistance(c(y1, 0), c(x2, a), type='GreatCircle')
			xy.wrong.lat <- pointDistance(c(0, x1), c(a, y2), type='GreatCircle')
			signch.lat.wrong.lon <- pointDistance(c(0, -y1), c(a, y2), type='GreatCircle')
			signch.lon.wrong.lat <- pointDistance(c(-x1, 0), c(x2, a), type='GreatCircle')
			
			signch.lat.xy.wrong.lon <- pointDistance(c(-y1, 0), c(x2, a), type='GreatCircle')
			signch.lon.xy.wrong.lat <- pointDistance(c(0, -x1), c(a, y2), type='GreatCircle')
			out <- c(difference, xy.exchange, signch.lat, signch.lon, signch.latlon, xy.signch.lat, xy.signch.lon, xy.signch.latlon, wrong.lon, wrong.lat, xy.wrong.lon, xy.wrong.lat, signch.lat.wrong.lon, signch.lon.wrong.lat, signch.lat.xy.wrong.lon, signch.lon.xy.wrong.lat)
			return(out)
		}
}

.errorGenerate <- function(pp){
		x1 <- pp[1]
		y1 <- pp[2]

		if(any(is.na(c(x1, y1)))){
			out <- NA
		} else {	
			no.change <- c(x1, y1)
			xy.exchange <- c(y1, x1)
			signch.lat <- c(x1, -y1)
			signch.lon <- c(-x1, y1)
			signch.latlon <- c(-x1, -y1)
			
			xy.signch.lat <- c(y1, -x1)
			xy.signch.lon <- c(-y1, x1)
			xy.signch.latlon <- c(-y1, -x1)
			
			out <- rbind(no.change, xy.exchange, signch.lat, signch.lon, signch.latlon, xy.signch.lat, xy.signch.lon, xy.signch.latlon)
			rownames(out) <- c("no.change", "xy.exchange", "signch.lat", "signch.lon", "signch.latlon", "xy.signch.lat", "xy.signch.lon", "xy.signch.latlon")
			return(out)
		}
}

	
		