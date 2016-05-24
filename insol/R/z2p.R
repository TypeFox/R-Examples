z2p <-
function(z,P0=101325,T0=288.15) {
	if (nargs() < 1 ) {cat("USAGE: z2p(z [,P0,T0]) \ units = m, Pa, K"); return()}
	Earth_G = 9.80665  		# Acceleration due to gravity (m s-2)
	EarthR = 6.3756766E6    # Average earths radius (m)
	Md = 28.966        		# Molecular weight of dry air
	R_star = 8.3145        	# Universal gas constant J/molK
	stlapse = -0.0065  		# standard lapse rate K/m
	H1 = (EarthR * z) /(EarthR + z)  
	HB = 0.0
	zp = P0*(T0/(T0+stlapse*(H1-HB)))**((Earth_G*Md)/(R_star*stlapse*1000))
	zp = zp/100.0  
	return(zp)
}

