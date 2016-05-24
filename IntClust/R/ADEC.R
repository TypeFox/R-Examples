ADEC<-function(List,distmeasure="tanimoto",normalize=FALSE,method=NULL,t=10,r=NULL,nrclusters=NULL,clust="agnes",linkage="ward",alpha=0.625,ResampleFeatures=TRUE){
	
	if(length(nrclusters)==1 & ResampleFeatures==TRUE){
		Clustering= ADECa(List=List,distmeasure=distmeasure,normalize=normalize,method=method,t=t,r=r,nrclusters=nrclusters,clust=clust,linkage=linkage,alpha=alpha)
	}
	else if(length(nrclusters)>1 & ResampleFeatures==FALSE){
		Clustering=ADECb(List=List,distmeasure=distmeasure,normalize=normalize,method=method,nrclusters=nrclusters,clust=clust,linkage=linkage,alpha=alpha)
	}
	else if(length(nrclusters)>1 & ResampleFeatures==TRUE){
		Clustering=ADECc(List=List,distmeasure=distmeasure,normalize=normalize,method=method,t=t,r=t,nrclusters=nrclusters,clust=clust,linkage=linkage,alpha=alpha)
		
	}
	else{
		stop("Cannot detect which version is requested. Please specify the desired number of clusters and whether ResampleFeatures is TRUE or FALSE.")
	}
	return(Clustering)
	
}