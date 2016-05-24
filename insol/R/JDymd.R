JDymd <-
function(year,month,day,hour=12,minute=0,sec=0) {
	if (nargs() < 3 ) {cat('USAGE: declination(year,month,day,hour=12,minute=0,sec=0) \n'); return()}
# 	valid 1901 to 2099
	hour=hour+minute/60+sec/60
	jd=367*year - (7*(year+(month+9)%/%12))%/%4 + (275*month)%/%9+day+1721013.5 + hour/24
return(jd)
}

