standardout<-function(Rdetector,iicc){
	
	Scores<-labels<-lapply(seq(1, length(Rdetector[[1]]), 1), function(x){lapply(seq(1, length(Rdetector), 1), function(y){
																								  vector(mode='numeric',length=length(Rdetector[[y]][[x]]))})})
	k<-length(iicc$DNA[[1]])
	col<-ncol(iicc$Transcriptionfactor) 
	
	for (jj in c(1:length(Rdetector))){
 	    for (i in c(1:length(Rdetector[[jj]]))) {
			Scores[[i]][[jj]]<-Rdetector[[jj]][[i]]
			for (p in c(1:length(iicc$position))){
				pos<-iicc$position[p]
				for (j in c(1:(nrow(iicc$Transcriptionfactor)-1))){
					r<-(k-ncol(iicc$Transcriptionfactor)+1)*(j-1)
					rr<-which(Scores[[i]][[jj]][r+c((pos-col):(pos+col-1))]==max(Scores[[i]][[jj]][r+c((pos-col):(pos+col-1))]))
					 labels[[i]][[jj]][r+ pos-col-1+rr]<-1
					 #labels[[i]][[jj]][r+ pos]<-1
					
				}
			}
		}
	}

	
y<-list(Scores=Scores,labels=labels)
	
	
return(y)

}
