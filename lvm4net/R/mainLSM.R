mainLSM<-function(psi2T, xiT, EZ, VZ, Y, xi, psi2, s2){
		
		SI4VZ <- solve(diag(nrow(VZ)) + 4 * VZ)
		
		cont <- dist(EZ %*% chol(SI4VZ))^2
		
		xiT <- upXiT(psi2T, xiT, SI4VZ, Y, cont, xi, psi2)
		psi2T <- upPsi2T(psi2T, xiT, SI4VZ, cont, psi2)

		VZ <- upSigmaT(psi2T, xiT, SI4VZ, Y, cont, EZ, s2)
		EZ <- upZT(psi2T, xiT, VZ, EZ, Y, s2)
			
return(list(xiT = xiT, psi2T = psi2T, lsmVZ = VZ, lsmEZ = EZ))

}