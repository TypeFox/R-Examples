sim.bd.taxa <-function(n,numbsim,lambda,mu,frac=1,complete=TRUE,stochsampling=FALSE){
	out<-lapply(1:numbsim,sim.bd.taxa.help,n=n,lambda=lambda,mu=mu,frac=frac,complete=complete,stochsampling=stochsampling)
	out1<-lapply(out,function(x){ x[[1]][[1]]})
	out2<-sapply(out,function(x){ x[[2]]})
	for (i in 1:length(out1)){
		out1[[i]]$root.edge<-out2[i]-max(getx(out1[[i]],sersampling=1)[,1])
	}
	#out3<-list(out1,out2)
	#out3
	out1
	}
