## Equation M.Groundspeed ##
NCEP.M.Groundspeed <- function(u, v, direction, p.airspeed,...){
		deg2rad = pi/180            # Conversion degrees to radians
		rad2deg = 180/pi		    # Conversion radians to degrees
	wdir <- ifelse(u==0 & v==0, direction, atan2(u,v)*rad2deg)	## Calculate wind direction
	wspd <- sqrt(u^2 + v^2)		## Calculate wind speed
	theta <- wdir-direction		## Calculate angle between wind and bird

	## Calculate tail and side wind components and the preferred groundspeed ##
	tailwind <- wspd*(cos(theta*deg2rad))
	sidewind <- wspd*(sin(theta*deg2rad))
	p.groundspeed <- p.airspeed + tailwind
	
	## Calculate forward and sideways movement ##
		side.move <- 0
		forward.move <- p.groundspeed
		
		## Calculate the bird's airspeed and groundspeed ##
		airspeed <- (sqrt(p.groundspeed^2 + wspd^2 - 2 * p.groundspeed * wspd * (cos(theta * deg2rad)))) 
		groundspeed <- p.groundspeed
		
		## Calculate the wind profit ##
		fa <- groundspeed - airspeed

	## Put all variables in a data.frame and return ##
	move <- data.frame(fa, forward.move, side.move, tailwind, sidewind, airspeed, groundspeed)
	return(move)
}