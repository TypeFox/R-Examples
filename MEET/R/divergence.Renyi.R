divergence.Renyi <-function(training.set,pmX,pmXY,q,correction){
	
	numpos<-1:ncol(training.set)					
	numfila<-1:4
	
	VarD<-((4-1)^2)/(2*((nrow(training.set))^2)*(log(2,base=exp(1)))^2)
	ErrorMI<-(2*VarD)^(1/2)
	
	
	r<-sapply(numpos,function(pos){
		sapply(numpos,function(position){
            #sha de revisar funció C
            
          rr<-matrix(0,4,4)
          out<-.C("powRenyi",
                            numfila=as.integer(4),
                            numpos=as.integer(ncol(training.set)),
                            pmX1=as.double(as.vector(pmX[,pos])),
                            pmX2=as.double(as.vector(pmX[,position])),
                            pmXY=as.double(as.vector(pmXY[[pos]][[position]])),
                            q=as.double(q),
                            rr=as.double(rr))
            
          z<-matrix(out$rr,4,4)
		  z[which(is.na(z))]<-0
		  sum(z)
		  })
	    })
    
	 rr<-(1/(q-1))*log2(r)
    
    # sample finite error correction
    
	
	if (q<1) {
		
		rr[rr+abs(ErrorMI)> slot(correction,"Herror")[nrow(training.set)]]<-slot(correction,"Herror")[nrow(training.set)]
		rr[rr-abs(ErrorMI)<0]<-0
	
	}else{
	 
		X<-sapply(numpos,function(pos){
				  sapply(numpos,function(position){
                        rr<-matrix(0,4,4)
                        out<-.C("powX2",
                                        numfila=as.integer(4),
                                        numpos=as.integer(ncol(training.set)),
                                        pmX1=as.double(as.vector(pmX[,pos])),
                                        pmX2=as.double(as.vector(pmX[,position])),
                                        pmXY=as.double(as.vector(pmXY[[pos]][[position]])),
                                        q=as.double(q),
                                        rr=as.double(rr))
                        z<-matrix(out$rr,4,4)
                        sum(z)
                        })
				  })
		 
        if ((length(which(rr+abs(ErrorMI)>(0.5*X))))!=0){
            
            rr[which(rr+abs(ErrorMI)>(0.5*X))]<-0.5*X[which(rr+abs(ErrorMI)>(0.5*X))]
        }
        
        if ((length(which(rr-abs(ErrorMI)<0)))!=0){
            
            rr[rr-abs(ErrorMI)<0]<-0
        }
    }
    return(rr)
}


