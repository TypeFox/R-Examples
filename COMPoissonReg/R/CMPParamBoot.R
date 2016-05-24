CMPParamBoot <-
function(x, poissonest, betahat, nuhat, n=1000){
	# Generate 1000 samples, using betahat and nuhat from full dataset
	Ystar <- matrix(0,nrow=nrow(x),ncol=n)
	for (i in 1:n){
	   Ystar[,i] <- makeCMPdata(x,betahat,nuhat)
	}
	
	# Take each of the 1000 sample results along with the x matrix, and run CMP regression on it to generate new betas and nu
	CMPresult <- matrix(0,nrow=n,ncol=length(betahat)+1)
	for (i in 1:n){
	   CMPresult[i,] <- ComputeBetasAndNuHat(as.matrix(x),Ystar[,i],poissonest,nuinit=1,max=100)$par 
	}
	
	return(list(Ystar=Ystar, CMPresult=CMPresult))
}

