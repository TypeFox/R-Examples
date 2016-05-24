divergence.Shannon <-function(training.set,H,HXY,correction){
	
	numpos<-1:ncol(training.set)	
	
	VarD<-((4-1)^2)/(2*((nrow(training.set))^2)*(log(2,base=exp(1)))^2)
	ErrorMI<-(1*VarD)^(1/2)

    Herror<-slot(correction,"Herror")[nrow(training.set)]
    Hxx<-mi<-matrix(0,ncol(training.set),ncol(training.set))

    out<-.C("divergenceShannon",numpos=as.integer(ncol(training.set)),H=as.double(H),Hxx=as.double(Hxx), HXY=as.double(HXY),mi=as.double(mi),Herror=as.double(Herror),ErrorMI=as.double(ErrorMI))			
	
        
    MI<-matrix(out$mi,ncol(training.set),ncol(training.set))
    
  	return(MI)

}

