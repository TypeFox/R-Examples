polyTermNorm <- function(degree, polynomType, alpha, beta){
	coeff = polyCoeff(degree,polynomType,alpha,beta)
	cN <- rep(0,degree+1)

    if (degree == 0){cN[1] = coeff[1]}
	else if (degree == 1){
		cN[1] = coeff[2]
		cN[2] = coeff[1]		
	} else {
		cNM1 = rep(0,degree)
		cNM2 = rep(0,degree-1)
		
		cNM1 <- polyTermNorm(degree-1,polynomType,alpha,beta)#,cNM1)
		cNM2 <- polyTermNorm(degree-2,polynomType,alpha,beta)#,cNM2)
		
		cN[1] = coeff[2]*cNM1[1] + coeff[3]*cNM2[1]
		for (ii in 0:(degree-1)){
			cN[ii+2]=coeff[1]*cNM1[ii+1]+coeff[2]*cNM1[ii+2]+coeff[3]*cNM2[ii+2]
		}
		cN[degree] = coeff[1]*cNM1[degree-1] + coeff[2]*cNM1[degree]
		cN[degree+1] = coeff[1]*cNM1[degree]
	}	
	return(cN)
}