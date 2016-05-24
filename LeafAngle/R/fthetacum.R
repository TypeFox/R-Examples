`fthetacum` <-
function(angles, degrees=TRUE, ...){

	ps <- c()
	for(i in seq_along(angles)){
		upangle <- ifelse(degrees, angles[i], angles[i]*pi/180)
		ps[i] <- integrate(ftheta,lower=0,upper=upangle, degrees=degrees, ...)[[1]]
	}
ps
}

