declination <-
function (jd){
	if (nargs() < 1 ) {cat('USAGE: declination(jd) \n'); return()}
	jdc=(jd - 2451545.0)/36525.0
	sec = 21.448 - jdc*(46.8150 + jdc*(0.00059 - jdc*(0.001813)))
	e0 = 23.0 + (26.0 + (sec/60.0))/60.0  
	oblcorr = e0 + 0.00256 * cos(radians(125.04 - 1934.136 * jdc))  
	l0 = 280.46646 + jdc * (36000.76983 + jdc*(0.0003032))
	l0 = (l0-360*(l0%/%360))%%360
	gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
	gmas=radians(gmas)
	seqcent = sin(gmas) * (1.914602 - jdc * (0.004817 + 0.000014 * jdc)) + 
		sin(2*gmas) * (0.019993 - 0.000101 * jdc) + sin(3*gmas) * 0.000289
	suntl = l0 + seqcent
	sal = suntl - 0.00569 - 0.00478 * sin(radians(125.04 - 1934.136 * jdc))
	delta = asin( sin(radians(oblcorr))*sin(radians(sal)) )
	return(degrees(delta))
}

