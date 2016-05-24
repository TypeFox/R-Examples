TrendWeights = function(evs.sum,parent.sum) {
	nvar = dim(evs.sum)[1]
	beta = rep(NA,nvar)
	for(i in c(1:nvar)){
		D0 = parent.sum[i,3]
		D1 = parent.sum[i,2]
		D2 = parent.sum[i,1]
		ND = D0+D1+D2
		
		H0 = evs.sum[i,3]
		H1 = evs.sum[i,2]
		H2 = evs.sum[i,1]
		NH = H0+H1+H2
		
		N0 = D0+H0
		N1 = D1+H1
		N2 = D2+H2
		
		N = N0+N1+N2
		
		if(N > N0) {
			beta[i] = (sqrt(N) * (N*(D1+2*D2) - ND*(N1+2*N2)))/sqrt(ND*NH*(N*(N1+4*N2) - (N1+2*N2)^2))
			#	beta[i] = (N*(N*(D1+2*D2) - ND*(N1+2*N2))^2) / (ND*NH*(N*(N1+4*N2) - (N1+2*N2)^2))
		}else{
			beta[i] = -1  ##weights for de novo mutations
		}
		
		
	}
	return(beta)
}	
