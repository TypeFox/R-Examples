#'@title Returns Schmidt number for a specific gas at a given temperature
#'@description 
#'Schmidt number is temperature dependant, and is the ratio of the kinematic viscosity of water 
#'to a diffusion coefficient. Coefficients are included for He, O2, CO2, CH4, SF6, N2O, Ar, and N2.\cr
#'
#'@usage
#'getSchmidt(temperature, gas)
#'
#'@param temperature Numeric vector of water temperatures in deg. Celsius
#'@param gas String for gas code. Valid inputs include: He, O2, CO2, CH4, SF6, N2O, Ar, and N2
#'@return Schmidt number (unitless)
#'@note Temperature range is only valid from 4-35 deg Celsius
#'@keywords methods math
#'@references
#'Raymond, Peter A., Christopher J. Zappa, David Butman, Thomas L. Bott, Jody Potter, Patrick Mulholland, 
#'Andrew E. Laursen, William H. McDowell, and Denis Newbold. \emph{Scaling the gas transfer velocity and hydraulic 
#'geometry in streams and small rivers}. Limnology & Oceanography: Fluids & Environments 2 (2012): 41-53.
#'@author
#'Jordan S. Read
#'@examples 
#'getSchmidt(temperature=12, gas="O2")
#'@export
getSchmidt	<-	function(temperature, gas){
	
	range.t	<-	c(4,35) # supported temperature range
	
	Schmidt	<-	data.frame(
		"He"=c(368,-16.75,0.374,-0.0036),
		"O2"=c(1568,-86.04,2.142,-0.0216),
		"CO2"=c(1742,-91.24,2.208,-0.0219),
		"CH4"=c(1824,-98.12,2.413,-0.0241),
		"SF6"=c(3255,-217.13,6.837,-0.0861),
		"N2O"=c(2105,-130.08,3.486,-0.0365),
		"Ar"=c(1799,-106.96,2.797,-0.0289),
		"N2"=c(1615,-92.15,2.349,-0.0240)
	)
		
	obsT <- is.finite(temperature) # logical for observed (not NA or NaN [or Inf or -Inf]) -RDB
		
	if (!is.character(gas)){stop(paste('gas must be a character. was given as',gas))}
	if (length(gas)>1){stop("only one gas can be specified for this version")}
	if (!any(names(Schmidt)==gas)){stop(paste(gas,'not found in list of coded gasses'))}
	if (any(temperature[obsT] < range.t[1] | temperature[obsT] > range.t[2])){ # This logical threw an error if any temperature were NA (or NaN, etc.) -RDB
		warning("temperature out of range")
	}
	A	<-	unlist(Schmidt[gas])[1]
	B	<-	unlist(Schmidt[gas])[2]
	C	<-	unlist(Schmidt[gas])[3]
	D	<-	unlist(Schmidt[gas])[4]

	Sc = as.numeric(A+B*temperature+C*temperature^2+D*temperature^3)

	return(Sc)
}