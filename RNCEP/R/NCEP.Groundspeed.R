## Equation Groundspeed ##
NCEP.Groundspeed <- function(u, v, direction, groundspeed,...){
		deg2rad = pi/180            # Conversion degrees to radians
		rad2deg = 180/pi		    # Conversion radians to degrees
	wdir <- ifelse(u==0 & v==0, direction, atan2(u,v)*rad2deg)	## Calculate wind direction
	wspd <- sqrt(u^2 + v^2)		## Calculate wind speed
	theta <- wdir-direction		## Calculate angle between wind and bird

	## Calculate the wind profit ##
	fa <- groundspeed - (sqrt(groundspeed^2 + wspd^2 - 2 * groundspeed * wspd * (cos(theta * deg2rad))))

	## Calculate the forward and side movements ##
	tailwind <- wspd*(cos(theta*deg2rad))
	sidewind <- wspd*(sin(theta*deg2rad))
	side.move <- 0
	forward.move <- groundspeed
	
	## Calculate the bird's airspeed and groundspeed ##
	airspeed <- (sqrt(groundspeed^2 + wspd^2 - 2 * groundspeed * wspd * (cos(theta * deg2rad))))
	groundspeed <- groundspeed
	
	## Put all variables in a data.frame and return ##
	move <- data.frame(fa, forward.move, side.move, tailwind, sidewind, airspeed, groundspeed)

return(move)
}