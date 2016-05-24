detector_2nOrdre_init <-function(training.set, val.set, iicc) {
	 	
	    iicc$D <- NULL; iicc$correction_1rOrdre <- NULL
	    iicc$correction_1rOrdre_m1 <- NULL;iicc$HX<-NULL
	    iicc$Mperfil<- NULL; iicc$HXmax<-NULL
    
        
		Prob 	<- as.numeric(iicc$background)
    
		q<-iicc$q		
		correction 	<-  correction.entropy(q,nrow(training.set), 1, iicc)
		correction_m1	<-  correction.entropy(q,nrow(training.set)+1 , 1, iicc)		
		Herror<-slot(correction,"Herror")
		HXmax<-slot(correction,"Herror")[nrow(training.set)]
	
        outprob<-.C("probabilitycouple",background=as.double(Prob),Probtrans=double(16))
        prob_parella<-outprob$Probtrans

		pmX     <- probability(training.set, Prob)
		pmXY    <- joint.probability(training.set,Prob,prob_parella)
		
		HXY	<- entropy.joint(pmXY,q,iicc)
	
		H<-switch(iicc$classentropy, "Shannon"=entropy.Shannon(pmX),"Renyi"=entropy.Renyi(pmX,q))

		Rcorregida<-redundancy(H,HXmax,Herror,finite=nrow(training.set),TRUE)	
	
		PE_perfil<-sqrt(as.vector(Rcorregida)%*%t(as.vector(Rcorregida)))
	
	    D<-switch(iicc$classentropy, "Shannon"=divergence.Shannon(training.set,H,HXY,correction),"Renyi"=divergence.Renyi(training.set,pmX,pmXY,q,correction))

		iicc 	<- c(iicc, q = q, probparella=list(prob_parella),D=list(D), correction_1rOrdre = correction, correction_1rOrdre_val = correction_m1,Mperfil=list(PE_perfil),Hmax=HXmax,HX=list(H) )
        return(iicc)
		
}

