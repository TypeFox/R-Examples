polyNorm <- function(degree,x,polynomType,alpha,beta){
	coeff <- polyCoeff(degree,polynomType,alpha,beta)
	if (degree == 0){
		value <- coeff[1]
	} else if (degree == 1){value <- (coeff[1]*x+coeff[2])*polyNorm(degree-1,x,polynomType,alpha,beta)
	} else {
		value <- (coeff[1]*x + coeff[2])*polyNorm(degree-1,x,polynomType,alpha,beta)
		value <- value + coeff[3]*polyNorm(degree-2,x,polynomType,alpha,beta)		
	}
	return (value)
}
