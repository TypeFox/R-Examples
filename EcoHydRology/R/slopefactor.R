slopefactor <- function(lat,Jday,slope,aspect){
# slopefactor: adjusts solar radiation for land slope and aspect relative to the sun, 1=level ground
#lat: latitdue [rad]
#Jday: Julian date or day of the year [day]
#slope: slope of the ground [rad]
#aspect: ground aspect [rad from north]
	SolAsp <- rep(pi, length(Jday))  # Average Solar aspect is binary - either north (0) or south (pi) for the day
	SolAsp[which(lat - declination(Jday) < 0)] <- 0   # 
	SF <- cos(slope) - sin(slope)*cos(aspect-(pi-SolAsp))/tan(solarangle(lat,Jday))
	SF[which(SF < 0)] <- 0  ## Slope factors less than zero are completely shaded

	return( SF )
}

