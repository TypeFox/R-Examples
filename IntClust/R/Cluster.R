Cluster<-function(Data,type=c("data","dist"),distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",alpha=0.625,gap=TRUE,maxK=50,StopRange=FALSE){	

	#STEP 1: Distance Matrices
	type<-match.arg(type)
	if(type=="data"){
		DistM=Distance(Data,distmeasure,normalize,method)
		if(StopRange==FALSE & !(0<=min(DistM) & max(DistM)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			DistM=Normalization(DistM,method="Range")
		}
	}
	else{
		DistM=Data
		if(StopRange==FALSE  & !(0<=min(DistM) & max(DistM)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			DistM=Normalization(DistM,method="Range")
		}
	}
	
	#STEP 2: Hierarchical Clustering with Ward Link
	
	Clust=cluster::agnes(DistM,diss=TRUE,method=linkage,par.method=alpha)
	
	func = function(x,k){
		return(list(cluster=stats::cutree(Clust,k=k)  )  )
	}
	
	#with optional gap statistic
	if(gap==TRUE){
		Clust_gap = cluster::clusGap(Data,FUNcluster=func,K.max=maxK,B=500)
		gapdata = as.data.frame(Clust_gap$Tab)
		
		k1 = cluster::maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"firstSEmax")
		k2 = cluster::maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"globalSEmax")
		k3 = cluster::maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"firstmax")
		k4 = cluster::maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"globalmax")
		k5 = cluster::maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"Tibs2001SEmax")
		
		k = data.frame(firstSEmax=k1,globalSEmax=k2,firstmax=k3,globalmax=k4,Tibs2001SEmax=k5)
		
		out = list(DistM=DistM,Clust=Clust,Clust_gap=Clust_gap,gapdata=gapdata,k=k)
	}
	else{
		out=list(DistM=DistM,Clust=Clust)
	}
	attr(out,'method')<-'Single Clustering'
	return(out)
}
