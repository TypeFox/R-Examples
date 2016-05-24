## Equation Tailwind ##
NCEP.Tailwind <- function(u, v, direction, airspeed=NA,...){
		deg2rad = pi/180            # Conversion degrees to radians
		rad2deg = 180/pi		    # Conversion radians to degrees
	wdir <- ifelse(u==0 & v==0, direction, atan2(u,v)*rad2deg)
	wspd <- sqrt(u^2 + v^2)
	theta <- wdir-direction		## Calculate angle between wind and bird
	## Calculate wind profit ##
	fa <- (wspd*cos((theta)*deg2rad))
	
	## Calculate the forward and side movement ##
	tailwind <- wspd*(cos(theta*deg2rad))
	sidewind <- wspd*(sin(theta*deg2rad))
	side.move <- sidewind
	forward.move <- airspeed + tailwind

	## Calculate the bird's ground and air speeds
	airspeed <- airspeed
	groundspeed <- sqrt(forward.move^2 + side.move^2) 
	
	## Put all variables in a data.frame and return ##
	move <- data.frame(fa, forward.move, side.move, tailwind, sidewind, airspeed, groundspeed)

return(move)
}