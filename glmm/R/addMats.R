addMats <-
function(matList){
		 if (! is.list(matList))
        matList <- list(matList)
        
        ncells<-length(as.vector(matList[[1]]))
        nr<-nrow(matList[[1]])
        matAsVec<-lapply(matList,as.vector)
        
        out<-rep(0,ncells)

        for(d in 1:ncells){
        	thing<-lapply(matAsVec,"[[",d)
        	thing<-unlist(thing)
        	out[d]<-sum(thing)
        }
	matrix(data=out,nrow=nr)
}

addVecs <-
function(vecs){
	 if (! is.list(vecs))
        vecs <- list(vecs)

	dvecs <-length(vecs[[1]])
	out<-rep(0,dvecs)
	
	
	for(d in 1:dvecs){
		thing<-lapply(vecs,"[[",d)
		thing<-unlist(thing)
		out[d]<-sum(thing)
	}
	out
}

catVecs <-
function(vecs){
	 if (! is.list(vecs))
        vecs <- list(vecs)

	lengths <-lapply(vecs,length)
	lengths<-unlist(lengths)
	out<-rep(0,sum(lengths))
	starthere<-1
	
	for(d in 1:length(vecs)){
		endhere<-starthere+lengths[d]-1
		out[c(starthere:endhere)]<-unlist(vecs[[d]])
		starthere<-starthere+lengths[d]
	}
	out
}
