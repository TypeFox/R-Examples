# Download geographic data and return as R object
# Author: Robert J. Hijmans, r.hijmans@gmail.com
# License GPL3
# Version 0.9
# October 2008

# compare some of this with geonames package, and perhaps use that instead. 



.alt <- function(lonlat) {
	lonlat <- .pointsToMatrix(lonlat)
	theurl <- paste("http://api.geonames.org/srtm3?lat=", lonlat[,2], "&lng=", lonlat[,2], "&username=demo", sep='')
	elevation <- scan(theurl, what='character', quiet=TRUE)
	if (elevation < -32000) { elevation <- NA }
	return(elevation)
}


.country <- function(lonlat, radius=0) {
	cnts <- ccodes()
	lonlat <- .pointsToMatrix(lonlat)

	res <- matrix(ncol=3,nrow=length(lonlat[,1]))
	for (i in 1:length(lonlat[,1])) {
		theurl <- paste("http://api.geonames.org/countryCode?lat=", lonlat[i,2], "&lng=", lonlat[i,1], "&radius=", radius, sep='')
		country <- scan(theurl, what='character', quiet=TRUE)
		if (length(country) > 1) { res[i,] <- c(NA,NA,NA)
		} else {
			rec <- subset(cnts, cnts[,3] == country) 
			if (length(rec) == 0) { res[i,] <- c(NA,NA,NA) 
			} else res[i,] <- rec
		}	
	}	
	colnames(res) <- c("NAME", "ISO3", "ISO2")
	return(res)
}


.adm <- function(lonlat, radius=0, maxrows=1) {
	lonlat <- .pointsToMatrix(lonlat)
	theurl <- paste("http://ws.geonames.org/countrySubdivision?lat=", lonlat[,1], "&lng=", lonlat[,2], "&radius=", radius, "&maxrows=", maxrows, sep='')
	subdivs <- scan(theurl, what='character', quiet=TRUE)
	return(subdivs)
}


#http://ws.geonames.org/findNearbyPlaceName?lat=47.3&lng=9 
#http://ws.geonames.org/findNearby?lat=47.3&lng=9 
#http://ws.geonames.org/findNearbyWikipedia?lat=47&lng=9


