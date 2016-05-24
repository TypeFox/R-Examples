SPU_E_diffcov <- function(gamma, n1, n2, cov1, cov2){
	if(gamma %% 2 == 1){
		e <- 0
	}
	if(gamma %% 2 == 0){
		ga.half <- gamma/2
		P <- 0
		for(d in seq(0, ga.half, 1)){
			P <- P + sum((diag(cov1))^d*(diag(cov2))^(ga.half - d))/(factorial(d)*factorial(ga.half - d)*n1^d*n2^(ga.half - d))
		}
		e <- factorial(gamma)/2^ga.half*P
	}
	return(e)
}