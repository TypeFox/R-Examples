gdist <-
function(lon.1, lat.1, lon.2, lat.2, units = 'nm', a = 6378137.0, b = 6356752.3142, verbose = FALSE)
{
#
# Calculate geodesic distance (in nm) between two points specified by latitude/longitude
# using Vincenty inverse formula for ellipsoids
# Reference: Direct and inverse solutions of geodesics on the ellipsoid with application
#  of nested equations.  Survey Review XXII, 176, April 1975.
#
# Inspired by: http://www.movable-type.co.uk/scripts/LatLongVincenty.html
#
#   DATE WRITTEN:  10 January 2005      LAST REVISED:   04 January 2010
#   AUTHOR:  John R. Wallace: Imap.for.R@gmail.com
#
#
        if(any(!is.finite(c(lon.1, lat.1, lon.2, lat.2))))
             return(NA)

	# a, b = major & minor semiaxes of the ellipsoid in meters
	# flat = flattening (a-b)/a
	# lat.1, lat.2 = geodetic latitude
	# L = difference in longitude
	
	
	rad <- pi/180
	lon.1 <- lon.1 * rad
	lat.1 <- lat.1 * rad
	lon.2 <- lon.2 * rad
	lat.2 <- lat.2 * rad
	
	flat <- (a - b)/a
	L <- lon.1 - lon.2
	U1 <- atan((1 - flat) * tan(lat.1))
	U2 <- atan((1 - flat) * tan(lat.2))
	
	lamda <- L
	lamda.old <- 2 * pi
	if(verbose)
		cat("\nStarting lamda =", lamda, "\n\n")
	
	i <- 1
	while(abs(lamda - lamda.old) > 1e-011) {
		sin.sigma <- sqrt((cos(U2) * sin(lamda))^2 + (cos(U1) * sin(U2) - sin(U1) * cos(U2) * cos(lamda))^2)
		cos.sigma <- sin(U1) * sin(U2) + cos(U1) * cos(U2) * cos(lamda)
		sigma <- atan2(sin.sigma, cos.sigma)
		sin.alpha <- (cos(U1) * cos(U2) * sin(lamda))/ifelse(sin(sigma) == 0, 1e-025, sin(sigma))
		cos2.alpha <- 1 - sin.alpha^2
		cos2.sigma.m <- cos(sigma) - (2 * sin(U1) * sin(U2))/ifelse(cos2.alpha == 0, 1e-025, cos2.alpha)
		C. <- (flat/16) * (cos2.alpha * (4 + flat * (4 - 3 * cos2.alpha)))
		if(verbose) {
			cat("sin.sigma =", sin.sigma, "\n")
			cat("cos.sigma =", cos.sigma, "\n")
			cat("sigma =", sigma, "\n")
			cat("sin.alpha =", sin.alpha, "\n")
			cat("cos2.alpha =", cos2.alpha, "\n")
			cat("cos2.sigma.m =", cos2.sigma.m, "\n")
			cat("C =", C., "\n")
			cat("lamda diff =", lamda - lamda.old, "\n")
		}
		lamda.old <- lamda
		lamda <- L + (1 - C.) * flat * sin.alpha * (sigma + C. * sin.sigma * (cos2.sigma.m + C. * cos.sigma * (-1 + 2 * cos2.sigma.m^2)))
		if(verbose)
			cat("New lamda =", lamda, "\n\n")
		if(i > 20) {
			warning("lamda did not converge")
			return(NA)
		}
		i <- i + 1
	}
	
	u2 <- (cos2.alpha * (a^2 - b^2))/b^2
	A <- 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
	B <- (u2/1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
	delta.sigma <- B * sin(sigma) * (cos2.sigma.m + (B/4) * (cos(sigma) * (-1 + 2 * cos2.sigma.m^2) - (B/6) * cos2.sigma.m * (-3 + 4 * sin(sigma)^2) * (
		-3 + 4 * cos2.sigma.m^2)))
	
	if(verbose) {
		alpha1 <- atan((cos(U2) * sin(lamda))/(cos(U1) * sin(U2) - sin(U1) * cos(U2) * cos(lamda)))
		cat("\nalpha1 =", alpha1/rad, "\n")
		alpha2 <- atan((cos(U1) * sin(lamda))/( - sin(U1) * cos(U2) + cos(U1) * sin(U2) * cos(lamda)))
		cat("alpha2 =", alpha2/rad, "\n\n")
		""
		cat("2*sigma.m =", acos(cos2.sigma.m), "\n")
		""
		cat("b =", b, "\n")
		cat("A =", A, "\n")
		cat("sigma (radians) =", sigma, "\n")
		cat("delta.sigma (radians) =", delta.sigma, "\n")
		cat("Distance: s = b * A * (sigma - delta.sigma)\n\n")
		""
	}
	
	s <- (b * A * (sigma - delta.sigma))

        switch(units, 
	      m = s,
             km = s/1000,
             nm = s/1852,
          miles = s/1609.344)



}


