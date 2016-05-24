JacobiPolyCoeff <- function(degree,alpha,beta){
	if (degree == 0){
		value <- c(1.0, 0.0, 0.0)
	} else if (degree == 1){
		value <- c((2.0+alpha+beta)/2.0, (alpha-beta)/2.0, 0.0)
	} else {
		a <- (2.0*degree-1.0+alpha+beta) * (2.0*degree+alpha+beta)   / (2.0*degree) / (degree+alpha+beta)
		b <- (alpha*alpha-beta*beta)     * (2*degree-1.0+alpha+beta) / (2.0*degree) / (2.0*degree-2.0+alpha+beta) / (degree+alpha+beta)
		c <- -(degree-1.0+alpha) * (degree-1.0+beta) * (2*degree+alpha+beta) / (degree) / (degree+alpha+beta) / (2.0*degree-2+alpha+beta)
		value <- c(a,b,c)	
	}
	return (value)
}