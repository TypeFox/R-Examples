## Equation Airspeed Movement ##
NCEP.Airspeed <- function(u, v, direction, airspeed,...){
		deg2rad = pi/180            # Conversion degrees to radians
		rad2deg = 180/pi		    # Conversion radians to degrees
	wdir <- ifelse(u==0 & v==0, direction, atan2(u,v)*rad2deg)	## Calculate wind direction
	wspd <- sqrt(u^2 + v^2)		## Calculate wind speed
	theta <- wdir-direction		## Calculate angle between wind and bird
	alpha <- pi - theta*deg2rad - asin((wspd*sin(theta*deg2rad))/airspeed)

## Calculate the wind profit ##
	fa <- wspd * cos(theta*deg2rad) + sqrt(airspeed^2 - (wspd * sin(theta*deg2rad))^2) - airspeed 


## Calculate the forward and side movements
	tailwind <- wspd*(cos(theta*deg2rad))
	sidewind <- wspd*(sin(theta*deg2rad))
	##
	side.move <- ifelse(is.na(fa), NA, 0)
	forward.move <- ifelse(is.na(fa), NA, sqrt(wspd^2 + airspeed^2 - 2*wspd*airspeed*cos(alpha)))

## Calculate the bird's ground and air speeds
	airspeed <- airspeed
	groundspeed <- forward.move

## Put all variables in a data.frame and return ##
	move <- data.frame(fa, forward.move, side.move, tailwind, sidewind, airspeed, groundspeed)
return(move)
}
