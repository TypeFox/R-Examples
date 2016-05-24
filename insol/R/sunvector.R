sunvector <-
function (jd,latitude,longitude,timezone) {
if (nargs() < 4 ) {cat("USAGE: sunvector(jd,latitude,longitude,timezone)\n values in jd, degrees, hours\n"); return()}
	omegar = hourangle(jd,longitude,timezone)
	deltar = radians(declination(jd))
	lambdar = radians(latitude)	
	svx = -sin(omegar)*cos(deltar)
	svy = sin(lambdar)*cos(omegar)*cos(deltar)-cos(lambdar)*sin(deltar)
	svz = cos(lambdar)*cos(omegar)*cos(deltar)+sin(lambdar)*sin(deltar) 
	return(cbind(svx,svy,svz))
}

