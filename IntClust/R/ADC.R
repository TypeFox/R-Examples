ADC<-function(List,distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",alpha=0.625){
	
	#Checking required data types and methods:
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	
	
	#Fuse variables into 1 data Matrix
	
	OrderNames=rownames(List[[1]])
	for(i in 1:length(List)){
		List[[i]]=List[[i]][OrderNames,]
	}
	
	AllData<-NULL
	for (i in 1:length(List)){
		if(i==1){
			AllData=List[[1]]
		}
		else{
			AllData=cbind(AllData,List[[i]])
		}
	}
	
	#Compute Distance Matrix on AllData
	
	AllDataDist=Distance(AllData,distmeasure,normalize,method)
	
	#Perform hierarchical clustering with ward link on distance matrix
	
	HClust = cluster::agnes(AllDataDist,diss=TRUE,method=linkage,par.method=alpha)		
	
	
	out=list(AllData=AllData,DistM=AllDataDist,Clust=HClust)
	attr(out,'method')<-'ADC'	
	return(out)
	
}

