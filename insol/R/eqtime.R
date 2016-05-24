eqtime <-
function(jd) {
	if (nargs() < 1 ) {cat("USAGE: eqtime(jd)\n"); return()}
	jdc=(jd - 2451545.0)/36525.0
	sec = 21.448 - jdc*(46.8150 + jdc*(0.00059 - jdc*(0.001813)))
	e0 = 23.0 + (26.0 + (sec/60.0))/60.0 
	ecc = 0.016708634 - jdc * (0.000042037 + 0.0000001267 * jdc)
	oblcorr = e0 + 0.00256 * cos(radians(125.04 - 1934.136 * jdc)) 
	y = (tan(radians(oblcorr)/2))^2
	l0 = 280.46646 + jdc * (36000.76983 + jdc*(0.0003032))
	l0 = (l0-360*(l0%/%360))%%360
	rl0 = radians(l0)
	gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
	gmas=radians(gmas)
	EqTime = y*sin(2*rl0)-2.0*ecc*sin(gmas)+4.0*ecc*y*sin(gmas)*cos(2*rl0)-
		0.5*y^2*sin(4*rl0)-1.25*ecc^2*sin(2*gmas)
	return(degrees(EqTime)*4)
}

