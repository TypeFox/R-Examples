simulat <-
function(x,pop=10000,nsamples=10,ssample=300)
{ 
{ 	
	ne<-round(x*pop)/100
	population<-rep(1:length(x),times=ne)
	population.t<-matrix(ncol=nsamples,nrow=length(population))
	for(i in 1:ncol(population.t))
		{
		population.t[,i]<-population
		}
	samples1<-apply(population.t,2,sample,ssample)
	samples<-apply(samples1,2,tabulate,nbins=length(x))
	rownames(samples)<-names(x)
	samples<-t(samples)
}
return(samples)
}

