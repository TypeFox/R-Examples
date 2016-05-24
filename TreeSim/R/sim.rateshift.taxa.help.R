sim.rateshift.taxa.help <-function(dummy,n,lambda,mu,frac,times,complete=TRUE,K,norm){
	out<-sim.rateshift.taxa.loop(n,1,lambda,mu,frac,times,complete,K,norm)
	out
	}