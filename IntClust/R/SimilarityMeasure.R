SimilarityMeasure<-function(List,nrclusters=NULL,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	
	if(class(List)!="list"){
		MatrixColors=List
	}	
	else{
		MatrixColors=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
	}
	
	#Compare every row to the first row
	Similarity=c()
	for(i in 1:dim(MatrixColors)[1]){
		Shared=0
		for(j in 1:dim(MatrixColors)[2])
			if(MatrixColors[i,j]==MatrixColors[1,j]){
				Shared=Shared+1
			}
		Similarity=c(Similarity,Shared/ncol(MatrixColors))
		
		
	}
	
	return(Similarity)
}
