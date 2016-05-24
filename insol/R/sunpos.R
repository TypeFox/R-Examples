sunpos <-
function(sunv) {
#  no refraction, center of disc
if (nargs() < 1 ) {cat("USAGE: sunpos(sunvector)\n 3D vector\n"); return()}
	azimuth = degrees(pi - atan2(sunv[,1],sunv[,2]  ) )
	zenith = degrees(acos(sunv[,3]))
	return(cbind(azimuth,zenith))
}

