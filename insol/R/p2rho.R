p2rho <-
function(Pz, TempK, RH){
	if (nargs() < 3 ) {cat("USAGE: result = function(Pressure, TempK, RH) \n"); return()}
	R_d = 287.04							# Gas constant for dry air J/KgK
	e_v_star = 0
	e_0 = 0
	P_d = 0
	rho_d = 0
	rho_v = 0
	e_v_star = wvapsat(TempK, 0)			# in hPa
	e_0 = e_v_star*RH/100.0 				# actual partial pressure of w. vapour,hPa
	P_d = Pz - e_0   						# partial pressure of dry air
	rho_d = 100.*P_d /(R_d*TempK) 			# Brutsaert 3.4 (*100  hPa to Pascal)
	rho_v = 0.622*100.*e_0/(R_d*TempK) 		# Brutsaert 3.5
	p2rho = rho_d + rho_v     				# density in kg/m^3
return(p2rho)
}

