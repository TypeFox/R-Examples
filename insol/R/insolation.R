insolation <-
function(zenith,jd,height,visibility,RH,tempK,O3,alphag) {
	if (nargs() < 8 ) { cat("USAGE: insolation(zenith,jd,height,visibility,RH,tempK,O3,alphag)"); return() } 
	if (min(tempK,na.rm=TRUE) < 153) { print("temperature should be in Kelvin"); return() }
	Isc = 1361.0   			# solar constant (Wm^(-2)) (1)
	zenith[zenith>90]=90
	theta = radians(zenith)
	ssctalb = 0.9  # single scattering albedo (aerosols)(Iqbal, 1983)
	Fc = 0.84      # ratio of forward to total energy scattered (Iqbal, 1983)
	Pz = z2p(height)
	Mr = 1.0/(cos(theta)+0.15*((93.885-zenith)^(-1.253)))
	Ma = Mr*Pz/1013.25
# 	#** Use Lowe(1977) Lowes polynomials for vapor pressure
	wvap_s =  wvapsat(tempK)
# 	#Wprec = 0.493*(RH/100.0)*wvap_s/tempK   #precipitable water in cm Leckner (1978)
	Wprec = 46.5*(RH/100.0)*wvap_s/tempK  #Prata 1996
	rho2 = sunr(jd)
	TauR = exp((-.09030*(Ma^0.84) )*(1.0+Ma-(Ma^1.01)) )
	TauO = 1.0-( ( 0.1611*(O3*Mr)*(1.0+139.48*(O3*Mr))^(-0.3035) ) - 
			0.002715*(O3*Mr)*( 1.0+0.044*(O3*Mr)+0.0003*(O3*Mr)^2 )^(-1))
	TauG = exp(-0.0127*(Ma^0.26))
	TauW = 1.0-2.4959*(Wprec*Mr)*( (1.0+79.034*(Wprec*Mr))^0.6828 + 
			6.385*(Wprec*Mr) )^(-1)
	TauA = ( 0.97-1.265*(visibility^(-0.66)) )^(Ma^0.9)   #Machler, 1983
	TauTotal = TauR*TauO*TauG*TauW*TauA   
	In = 0.9751*(1/rho2)*Isc*TauTotal
	tauaa = 1.0-(1.0-ssctalb)*(1.0-Ma+Ma^1.06)*(1.0-TauA)
	Idr = 0.79*(1/rho2)*Isc*cos(theta)*TauO*TauG*TauW*tauaa*0.5*(1.0-TauR)/(1.0-Ma+Ma^(1.02))
	tauas = (TauA)/tauaa
	Ida = 0.79*(1/rho2)*Isc*cos(theta)*TauO*TauG*TauW*tauaa*Fc*(1.0-tauas)/(1.0-Ma+Ma^1.02)
	alpha_atmos = 0.0685+(1.0-Fc)*(1.0-tauas)
	Idm = (In*cos(theta)+Idr+Ida)*alphag*alpha_atmos/(1.0-alphag*alpha_atmos)
	Id = Idr+Ida+Idm
	return(cbind(In,Id))
}

