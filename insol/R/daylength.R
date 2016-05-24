daylength <-
function (lat, long, jd, tmz){
	if (nargs() < 4 ) {cat("USAGE: daylength(latitude, longitude, jd, timezone) \n values in degrees, julian days, hours \n"); return()}
		EqTime = eqtime(jd)
		delta = declination(jd)
		tanlatdel = -tan(radians(lat)) * tan(radians(delta))
		tanlatdel[tanlatdel>1]=1
		omega = acos(tanlatdel)
		daylen = (2*omega)/(2*pi/24)
		stndmeridian = tmz*15
		deltaLatTime=long-stndmeridian
		deltaLatTime = deltaLatTime * 24/360 
		sunrise = 12*(1-omega/pi)-deltaLatTime-EqTime/60 
		sunset = 12*(1+omega/pi)-deltaLatTime-EqTime/60
		sunrise[omega==0]=NA
		sunset[omega==0]=NA
		return(round(cbind(sunrise,sunset,daylen),2))
}

