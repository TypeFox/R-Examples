sim.bd.taxa.help <-function(dummy,n,lambda,mu,frac=1,complete=TRUE,stochsampling=FALSE){
	out<-sim.bd.taxa.loop(n,1,lambda,mu,frac,complete,stochsampling)
	out
	}