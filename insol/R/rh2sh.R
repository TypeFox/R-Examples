rh2sh <-
function(RH, tempk, Pz, ice) {
	if (nargs() < 1 ) {cat("USAGE: result = rh2sh(RH,TemK,Pz[,ice=1])\n"); return()}
    if (min(tempk,na.rm=TRUE) < 153.0){ cat("temperature should be in Kelvin\n"); return()}
    if (missing(ice)) {ice = 0}
    R_d = 287.04						# Gas constant for dry air J/KgK
    R_v = 461.40						# Gas constant for water vapor J/KgK
	e_v_star = wvapsat(tempk,0) 
	if (ice > 0) {e_v_star = wvapsat(tempk,1) }
	e_0 = e_v_star*RH/100.0 			# actual partial pressure of w. vapour,hPa
	P_d = Pz - e_0   					# partial pressure of dry air
# 	R_Rdv = R_d/R_v						# Ratio Rd/Rv = 0.622
# 	q = (R_Rdv*e_0)/(p_d+R_Rdv*e_0) 	# specific humidity. Jacobson99 (2.27)
	rho_d = 100.*P_d /(R_d*tempk)   	# Brutsaert 3.4 (*100  hPa to Pascal)
	rho_v = 0.622*100.*e_0/(R_d*tempk)   					# Brutsaert 3.5
	rho = rho_d + rho_v       			# density in kg/m^3
	q = rho_v/rho        				# specific humidity air, Brutsaert 3.2
	q = q*1000.          				# q is in g/kg
return(q)
}

