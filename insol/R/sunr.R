sunr <-
function(jd) {
	jdc=(jd - 2451545.0)/36525.0
	ecc=0.016708634-jdc*(0.000042037+0.0000001267*jdc)
	gmas=357.52911+jdc*(35999.05029-0.0001537*jdc)
	gmasr=radians(gmas)
	seqc=sin(gmasr)*(1.914602-jdc*(0.004817+0.000014*jdc))+sin(2*gmas)*(0.019993-0.000101*jdc)+
		sin(3*gmasr)*0.000289
	sta=gmas+seqc
	sunrv=(1.000001018 * (1 - ecc^2)) / (1 + ecc * cos(radians(sta)))
	return(sunrv)
}

