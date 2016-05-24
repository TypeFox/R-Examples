Ell<-function(psi2T, xiT, SigmaT, ZT, Y){
	
	SI4SigmaT <- solve(diag(nrow(SigmaT)) + 4 * SigmaT)
	
	A <- log(1 + sqrt(det(SI4SigmaT)) * exp(xiT + psi2T / 2) * exp(- dist(ZT %*% chol(SI4SigmaT))^2))
	
	sum((xiT - 2 * sum(diag(SigmaT)) - as.matrix(dist(ZT)^2)) * Y) - 2 * sum(A) 

	}