wvapsat <-
function(tempk,ice) {
	if (nargs() < 1 ) {cat("USAGE: wvapsat(tempk [,ice])\n K, [0,1]"); return()}
	if (min(tempk,na.rm=TRUE) < 153.0){ print("temperature should be in Kelvin"); return()}
	if (missing(ice)) {ice = 0}
	tempcl = tempk
	a0 = 6984.505294	
	a1 = -188.9039310
	a2 = 2.133357675	
	a3 = -1.288580973e-2	
	a4 = 4.393587233e-5
	a5 = -8.023923082e-8
	a6 = 6.136820929e-11
	if (ice > 0) { 
		tempcl = tempk - 273.15
		a0 = 6.109177956
		a1 = 5.03469897e-1
		a2 = 1.886013408e-2
		a3 = 4.176223716e-4
		a4 = 5.824720280e-6
		a5 = 4.838803174e-8
		a6 = 1.838826904e-10
	}
	return( a0+tempcl*(a1+tempcl*(a2+tempcl*(a3+tempcl*(a4+tempcl*(a5+tempcl*a6)))))  )
}

