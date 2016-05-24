ExpProcess<-function(ExpMatrix,pert){

pertlist<-c()
sampleNum<-dim(ExpMatrix)[1]
    for (i in (1:dim(ExpMatrix)[1]))
      {
	    Numsample<-length(which(as.numeric(ExpMatrix[i,]) ==0))
		pertlist<-c(pertlist,(Numsample/(dim(ExpMatrix)[2])))
	  }
	 NewMatri<-cbind(pertlist,ExpMatrix)
	 NewMatri<- NewMatri[which(NewMatri[,1]<(1-pert)),2:(dim(NewMatri)[2])]
     
	return(NewMatri)
}
