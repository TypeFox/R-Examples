sim.bd.taxa.loop <-function(n,numbsim,lambda,mu,frac=1,complete=TRUE,stochsampling=FALSE){
	if (complete == TRUE) {
		phy <- sim2.bd.reverse(round(n/frac),numbsim,lambda,mu)
	} else if (stochsampling==FALSE) {
		phy2 <- sim2.bd.reverse(round(n/frac),numbsim,lambda,mu)
		phy1 <- reconstructed.taxa(phy2[[1]],(round(n/frac)-n))
		phy<-list(phy1,phy2[[2]])
	} else {
		phy <- sim2.bd.fast(n,numbsim,lambda,mu,frac)
	}
	phy
	}
