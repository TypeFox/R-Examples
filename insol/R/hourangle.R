hourangle <-
function(jd,longitude,timezone) {
	if (nargs() < 3 ) {cat("USAGE: hourangle(jd,longitude,timezone)\n julian day, degrees, hours. Return radians \n"); return()}
	hour = ((jd-floor(jd))*24+12) %% 24
	eqtime = eqtime(jd)
	stndmeridian = timezone*15	    			
	deltalontime = longitude-stndmeridian 		
	deltalontime = deltalontime * 24.0/360.0    
	omegar = pi*( ( (hour + deltalontime + eqtime/60)/12.0 ) - 1.0) 
	return(omegar)
}

