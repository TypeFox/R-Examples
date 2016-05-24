polyCoeff <- function(degree,polynomType,alpha,beta) {
	coeff <- c(0.0, 0.0, 0.0)
	if (polynomType == "HERMITE"){
		coeff <- HermitePolyCoeff(degree)
	} else if (polynomType == "JACOBI") {
	  coeff <- JacobiPolyCoeff(degree,alpha,beta)	  
	} else if (polynomType == "LEGENDRE"){ 
		coeff <- JacobiPolyCoeff(degree,0,0)
	} else if (polynomType == "LAGUERRE"){
		coeff <- LaguerrePolyCoeff(degree,alpha)	
	}
	return(coeff)
}