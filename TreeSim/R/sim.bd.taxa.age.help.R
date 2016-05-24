sim.bd.taxa.age.help <-function(dummy,n,lambda,mu,frac=1,age,mrca){
	out<-sim.bd.taxa.age.loop(n,1,lambda,mu,frac,age,mrca)[[1]]
	if (mrca == FALSE) {out$root.edge<-age-max(branching.times(out))}
	out
	}


