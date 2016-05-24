upPsi2T<-function(psi2T, xiT, SI4SigmaT, cont, psi2){
	
	1 / (1 / psi2 + 2 * sum(1 / (1 + exp(-psi2T / 2) / (sqrt(det(SI4SigmaT)) * exp(xiT) * exp(-cont)))))

	}