SNFa=function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,NN=20,mu=0.5,T=20,clust="agnes",linkage="ward",alpha=0.625,StopRange=FALSE){

	#Checking required data types and methods:
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	
	if(mu<0.3 | mu >0.8){
		message("Warning: mu is recommended to be between 0.3 and 0.8 for the SNF method. Default is 0.5.")
	}
	
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}
	
	
	#STEP 1: Distance Matrices
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		DistM=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize,method))
		DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		DistM=List
		DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
	}
	else{
		DistM=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
		
		OrderNames=rownames(DistM[[1]])
		for(i in 1:length(DistM)){
			DistM[[i]]=DistM[[i]][OrderNames,OrderNames]
		}
	}
	

	#STEP 2: Affinity Matrices
	
	AffM=lapply(seq(length(List)), function(i) SNFtool::affinityMatrix(DistM[[i]], NN, mu))
	
	#STEP 3: Fuse Networks Into 1 Single Network
	
	SNF_FusedM=SNFtool::SNF(AffM, NN, T)
	rownames(SNF_FusedM)=rownames(List[[1]])
	colnames(SNF_FusedM)=rownames(List[[1]])
	Dist=1-SNF_FusedM
	
	#STEP 4: Perform Hierarchical Clustering with WARD Link

	HClust = cluster::agnes(Dist,diss=TRUE,method=linkage,par.method=alpha)		
	
	
	#Output= list with the fused matrix and the performed clustering
	out=list(SNF_FusedM=SNF_FusedM,DistM=Dist,Clust=HClust)
	attr(out,'method')<-'SNF'
	return(out)
}
