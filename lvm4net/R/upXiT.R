upXiT<-function(psi2T, xiT, SI4SigmaT, Y, cont, xi, psi2){
	
	A <- sqrt(det(SI4SigmaT)) * exp(-cont) * exp(psi2T / 2)
	
	(xi+ psi2 * (sum(Y) - 2 * sum(1 / (1 + exp(-xiT) / A)) + xiT * 2 * sum(1 / (2 + exp(-xiT) / A+ exp(xiT) * A)))) / (1 + psi2 * 2 * sum(1 / (2 + exp(-xiT) / A + exp(xiT) * A)))
	
	}