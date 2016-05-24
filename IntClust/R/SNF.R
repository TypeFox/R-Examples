SNF<-function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,NN=20,mu=0.5,T=20,clust="agnes",linkage="ward",alpha=0.625,StopRange=FALSE,Version="SNFa"){
	if(Version=="SNFa"){
		Clustering=SNFa(List,type=type,distmeasure=distmeasure,normalize=normalize,method=method,NN=NN,mu=mu,T=T,clust=clust,linkage=linkage,alpha=alpha,StopRange=StopRange)
			
	}
	if(Version=="SNFb"){
		Clustering=SNFb(List,type=type,distmeasure=distmeasure,normalize=normalize,method=method,NN=NN,mu=mu,T=T,clust=clust,linkage=linkage,alpha=alpha,StopRange=StopRange)
		
	}
	if(Version=="SNFc"){
		Clustering=SNFc(List,type=type,distmeasure=distmeasure,normalize=normalize,method=method,NN=NN,mu=mu,T=T,clust=clust,linkage=linkage,alpha=alpha,StopRange=StopRange)
		
	}
	
	return(Clustering)
	
}	