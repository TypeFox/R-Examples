sim.bd.age.help <-function(dummy,age,lambda,mu,frac=1,mrca=FALSE,complete=TRUE,K){
	out<-sim.bd.age.loop(age,1,lambda,mu,frac,mrca,complete,K)[[1]]
	if (class(out)=="phylo" && mrca == FALSE) {
		out$root.edge<-age-max(getx(out,sersampling=1)[,1])}
	out
	}
