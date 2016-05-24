sim.bd.taxa.age.loop <-
function(n,numbsim,lambda,mu,frac=1,age,mrca=FALSE){
	phy<-sim2.bd.fast.age(n,numbsim,lambda,mu,frac,age,mrca)
	phy
}

