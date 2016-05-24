upZT<-function(psi2T, xiT, SigmaT, ZT, Y, s2){
			
	SI4SigmaT <- solve(diag(nrow(SigmaT)) + 4 * SigmaT)
	
	for (n in 1:nrow(ZT)){
		
		dnj <- t(ZT[n,] - t(as.matrix(ZT[-n,])))
		
		A <- 1 / sqrt(det(SI4SigmaT)) * exp(- xiT - psi2T / 2) * exp(apply(dnj, 1, function(x) t(x %*% SI4SigmaT %*% x)))
		
		f1znTo <- - 2 * SI4SigmaT %*% colSums(dnj / (1+A))
		f2znTo <- - 2 * sum(1 / (1 + A)) * SI4SigmaT + 
			 4 * SI4SigmaT %*% (t(dnj / (2 + 1 / A + A)) %*% dnj) %*% SI4SigmaT
		
		numZnT <- colSums(as.matrix(ZT[-n,]) * (Y[n,-n]+Y[-n,n])) - t(f1znTo) + t(ZT[n, ]) %*% f2znTo 
		denZnT <- (sum(Y[n,-n] + Y[-n,n]) + 1 / (2 * s2)) * diag(nrow(SigmaT)) + f2znTo
		
		ZT[n, ]<- numZnT %*% solve(denZnT)
	}

	ZT
	}