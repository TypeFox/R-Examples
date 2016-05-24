SatVaporDensity <- function(T_C){
	#	T_C	= Temperature [C]
	VP <- SatVaporPressure(T_C)
	return(round(VP/(0.462 * (T_C+273.15)), 4))
}

