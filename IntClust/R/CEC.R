CEC<-function(List,distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,t=10,r=NULL,nrclusters=NULL,weight=NULL,clust="agnes",linkage=c("flexible","flexible"),alpha=0.625,WeightClust=0.5,StopRange=FALSE,ResampleFeatures=TRUE){
	
	if(length(nrclusters)==1 & ResampleFeatures==TRUE){
		Clustering= CECa(List=List,distmeasure=distmeasure,normalize=normalize,method=method,t=t,r=t,nrclusters=nrclusters,weight=weight,clust=clust,linkage=linkage,alpha=alpha,WeightClust=WeightClust,StopRange=StopRange)
			
	}
	else if(length(nrclusters)>1 & ResampleFeatures==FALSE){
		Clustering=CECb(List=List,distmeasure=distmeasure,normalize=normalize,method=method,nrclusters=nrclusters,weight=weight,clust=clust,linkage=linkage,alpha=alpha,WeightClust=WeightClust,StopRange=StopRange)
	}		
	else if(length(nrclusters)>1 & ResampleFeatures==TRUE){
		Clustering=CECc(List=List,distmeasure=distmeasure,normalize=normalize,method=method,t=t,r=t,nrclusters=nrclusters,weight=weight,clust=clust,linkage=linkage,alpha=alpha,WeightClust=WeightClust,StopRange=StopRange)
			
	}
	else{    
		stop("Cannot detect which version is requested. Please specify the desired number of clusters and whether ResampleFeatures is TRUE or FALSE.")
	}

	return(Clustering)
}
	