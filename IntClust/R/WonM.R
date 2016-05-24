WonM=function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,nrclusters=seq(5,25,1),clust="agnes",linkage=c("flexible","flexible"),alpha=0.625,StopRange=FALSE){
	
	type<-match.arg(type)
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}
	
	
	#Step 1: Distance Matrices
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize,method))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
	}
	else{
		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		OrderNames=rownames(Dist[[1]])
		for(i in 1:length(Dist)){
			Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
		}
	}
	
	
	#Step 2: perform hierarchical clustering on both distance matrices

	HClustering=lapply(seq(length(List)),function(i) cluster::agnes(Dist[[i]],diss=TRUE,method=linkage[i],par.method=alpha))

	
	#Step 3: cut the dendrograms into a range of K values
	
	#Give 0 to pair belonging together, give 1 to a pair not belonging together : ==> Distances created otherwise similarities.
	ClusterMembers<-function(HClust,nrclusters){
		Temp=lapply(seq(length(nrclusters)),function(i) stats::cutree(HClust,nrclusters[i]))		
		CM=lapply(seq(length(nrclusters)),function(i) matrix(1,dim(List[[1]])[1],dim(List[[1]])[1]))
		
		clusters<-function(temp,cm){
			for(l in 1:length(temp)){
				label=temp[l]
				sameclust=which(temp==label)
				cm[l,sameclust]=0			
			}
			return(cm)
		}
		
		CM2=lapply(seq(length(nrclusters)),function(i) clusters(temp=Temp[[i]],cm=CM[[i]]))
		Consensus2=Reduce("+",CM2)
		return(Consensus2)
		
		
	}
	
	Consensus=lapply(seq(length(List)), function(i) ClusterMembers(HClustering[[i]],nrclusters))
	
	OverallConsensus=Reduce("+",Consensus)	
	OverallConsensus=as.matrix(OverallConsensus)
	rownames(OverallConsensus)=rownames(Dist[[1]])
	colnames(OverallConsensus)=rownames(Dist[[1]])
	OverallClustering=cluster::agnes(OverallConsensus,diss=TRUE,method="ward")
	
	out=list("Single Distances"=Dist,ClustSep=HClustering,DistM=OverallConsensus,Clust=OverallClustering)
	attr(out,'method')<-'WonM'
	return(out)
}
