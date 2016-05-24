correction.entropy <-function(q,p,long,iicc){
	
	Prob<-as.numeric(iicc$background)
	Herror<-sderror<-VAR<-vector("numeric",p)
	
    out<-.C("correctionentropy",
        nrowTFBS=as.double(p),
        q=as.double(q),
        Prob=as.double(Prob),
        classentropy=as.character(iicc$classentropy),
        Herror=double(p),
        VAR=double(p))
    
    VAR<-out$VAR
    Herror<- out$Herror
    sderror<-(long*VAR)^(1/2)
    resultat<-new("correction")
    slot(resultat,"Herror")<-Herror
    slot(resultat,"sderror")<-sderror
    return(resultat)
}

