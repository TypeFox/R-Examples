polyWeight <- function(degree,xi,ii,polynomType,alpha,beta){
	if (polynomType == "HERMITE"){
		value <- 2^(degree-1)*factorial(degree)*sqrt(pi)/degree^2		
		#print(value)
		a <- polyNorm(degree-1,xi[ii],polynomType,alpha,beta)
		value = value/a^2
		#print(a)
		#warning("Hermite")
	} else if ((polynomType == "JACOBI") | (polynomType == "LEGENDRE")) {
		constant = 1.0
		for (jj in 1:length(xi)){
			if (jj != ii){
				constant = constant*(xi[ii]-xi[jj])				
			}
		}
        a <- polyNorm(degree+1,xi[ii],polynomType,alpha,beta)
		value = (2.0*degree+2.0+alpha+beta)*gamma(degree+alpha+1.0)*gamma(degree+beta+1.0)*2.0^(alpha+beta)
        value = value/(degree+alpha+beta+1.0)/gamma(degree+alpha+beta+1.0)/factorial(degree+1)
        value = -value/constant/a 
	} else if (polynomType == "LAGUERRE") {
		a <- polyNorm(degree+1,xi[ii],polynomType,alpha,beta)
		value = xi[ii]/(degree+1.0)^2/a^2
	}
	return (value)
}