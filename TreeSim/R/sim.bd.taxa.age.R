sim.bd.taxa.age <-function(n,numbsim,lambda,mu,frac=1,age,mrca=FALSE){
	out<-lapply(1:numbsim,sim.bd.taxa.age.help,n=n,lambda=lambda,mu=mu,frac=frac,age=age,mrca=mrca)
	out
	}