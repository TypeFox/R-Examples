simulat.t <-
function(x,pop=1000,nsamples=10,ssample=300,
	percenta=TRUE,last)
{
{
	samples.t<-matrix(nrow=nrow(x)*nsamples,ncol=ncol(x)+1)
 	spp<-colnames(x)
 	colnames(samples.t)<-c("Sample",spp)
 	samples.t[,1]<-rep(1:nrow(x),each=nsamples)
 	a<-matrix(ncol=2,nrow=nrow(x))
 	for(i in 1:nrow(a))
 		{
 		a[i,1]<-i*nsamples-nsamples+1
 		a[i,2]<-i*nsamples
 		samples.t[c(a[i,1]:a[i,2]),-1]<-simulat(x[i,],pop=pop,
		nsamples=nsamples,ssample=ssample)
		}
	if(percenta==TRUE){
		percentages1<-percenta(samples.t[,-1],first=2,last=last)
		percentages<-matrix(nrow=nrow(x),ncol=ncol(x)*nsamples)
		colnames(percentages)<-rep(colnames(x),each=nsamples)
		for(i in 1:nrow(a))
			{
			percentages[i,]<-as.vector(percentages1[a[i,
				1]:a[i,2],])
			}
		}
}
if(percenta==TRUE){
	results<-list(samples.t,percentages)
	names(results)<-c("samples.t","percentages")
	return(results)
	}
else{
	return(samples.t)
	}
}

