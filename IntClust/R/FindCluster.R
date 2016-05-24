FindCluster<-function(List,nrclusters=NULL,select=c(1,1),fusionsLog=TRUE, WeightClust=TRUE,names=NULL){
	if(length(List)==1 & attributes(List[[1]])$method == "Weighted" & WeightClust==TRUE){
		T=List[[1]]$Clust
		attr(T,"method")="Single Clustering"
		List=list(T)
	}
	Matrix=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
	methodnr=select[1]
	clusternr=select[2]
	Comps=names(which(Matrix[methodnr,]==clusternr))
	return(Comps)
}	