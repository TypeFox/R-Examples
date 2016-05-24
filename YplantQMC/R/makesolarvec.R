	makesolarvec <- function(altitude, azimuth){
		k <- pi / 180
		azimuth <- azimuth - 180 # !!! Watch out; different azimuth definition.
		xdir <- sin(k * azimuth) * cos(k*altitude)
		zdir <- cos(k * azimuth) * cos(k*altitude)
		ydir <- -sin(k * altitude)
		if(length(azimuth) > 1)stop("Only support one light source, for now.")
		solarvec <- round(c(xdir,ydir,zdir,1),4)
	return(solarvec)
	}
