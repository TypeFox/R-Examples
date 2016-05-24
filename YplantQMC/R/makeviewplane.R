	# Returns vectors of x and y axes (list with x=c(x,y,z) and y=c(x,y,z))
	# of the 'viewplane'.
	# Taken from Yplant, procedure with the same name.
	makeviewplane <- function(azimuth, altitude){
		k <- pi/180
		x <- c(-cos(azimuth*k),sin(azimuth*k),0)
		y <- c(-sin(azimuth*k)*sin(altitude*k),
			   -cos(azimuth*k)*sin(altitude*k),
			   cos(altitude*k) )
	    z <- xprod(x,y)
	return(list(x=x,y=y,z=z))
	}