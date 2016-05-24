## Adriaan's Equation for Partial Compensation based on a proportion of the side wind speed ##
NCEP.PartialSpeed <- function(u, v, direction, airspeed, f=0.5,...){
	deg2rad = pi/180		 # Conversion degrees to radians
    rad2deg = 180/pi         # Conversion radians to degrees
    wdir <- ifelse(u==0 & v==0, direction, atan2(u,v)*rad2deg)
	wspd <- sqrt(u^2 + v^2)
    theta <- wdir - direction

## Calculate tail and side wind components ##
	tailwind <- wspd*(cos(theta*deg2rad))
	sidewind <- wspd*(sin(theta*deg2rad))
	

## Calculate the forward and side movements ##
	forward.move <- tailwind + sqrt(airspeed^2 - (f * sidewind)^2)
	side.move <- ifelse(is.na(forward.move), NaN, (1-f) * sidewind)

	## Calculate ground and air speeds ##
	airspeed <- airspeed
	groundspeed <- sqrt(forward.move^2 + ((1-f) * sidewind)^2)

## Calculate wind profit ##
	fa <- forward.move - airspeed 

## Put all variables in a data.frame and return ##
	move <- data.frame(fa, forward.move, side.move, tailwind, sidewind, airspeed, groundspeed)
return(move)
	}
