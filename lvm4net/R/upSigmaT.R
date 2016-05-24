upSigmaT<-function(psi2T, xiT, SI4SigmaT, Y, cont, ZT, s2){
	
	D<-nrow(SI4SigmaT)
	N<-nrow(Y)

	A <- 1 + 1 / sqrt(det(SI4SigmaT)) * exp(- xiT - psi2T / 2) * exp(cont)
	
	cbZT <- combn(nrow(ZT), 2)
	dij <- as.matrix(ZT[cbZT[1, ], ] - ZT[cbZT[2, ], ])
	
	f1s2To <- 8 * SI4SigmaT %*% (t(dij / A) %*% dij) %*% SI4SigmaT - 4 * sum(1 / A) * SI4SigmaT
	
	solve((1 / s2 + 4 * sum(Y) / N) * diag(D) + 2 * f1s2To / N)
	
}	