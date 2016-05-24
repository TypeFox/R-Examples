joint.probability <-function(training.set,Prob,Probtrans){
	
	yyy<-lapply(c(1:ncol(training.set)), function(i){
				 lapply(c(1:ncol(training.set)), function(ii){
						
						
                       outprob<- .C("jointprobability",
                                    trainingset=as.character(paste(training.set[,i],training.set[,ii],sep='')),
                                    indexA=as.integer(i),
                                    indexB=as.integer(ii),
                                    nrowTFBS=as.integer(nrow(training.set)),
                                    Probtrans=as.double(Probtrans),
                                    background=as.double(Prob),
                                    Prob=double(16))
                            
						yy<-matrix( data = outprob$Prob, nrow = 4, ncol = 4, byrow = TRUE,
                                    dimnames=list(c("A","T","C","G"),c("A","T","C","G")))

						})
				 })
	return(yyy)
	
}

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	



