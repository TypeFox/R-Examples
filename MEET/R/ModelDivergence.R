ModelDivergence<-function(iicc){
	
    require("MEET")	
	Prob<-iicc$background
	Prob<-as.numeric(Prob)
	training.set<-Factortrans<-matriu<-iicc$Transcriptionfactor
	#iicc$missing.fun=TRUE
	
	iicc$D <- iicc$correction_1rOrdre <- iicc$probparella<-NULL
	iicc$correction_1rOrdre_m1<-iicc$Mperfil<-iicc$q<- NULL
	iicc$HXmax<-iicc$HX<-iicc$q<-NULL
	iicc$Divergence<-iicc$interA<-iicc$interB<-iicc$classentropy<-NULL
	#iicc$missing.fun=TRUE
	q<-iicc$q<-iicc$parametersIdeal	
	if (q==1) {iicc$classentropy<-"Shannon"
		}else{
		iicc$classentropy<-"Renyi"  
	}
	
	iicc <- detector_2nOrdre_init(training.set, val.set=NULL, iicc)
    
	memory<-MImemory(iicc,Factortrans)
	iicc<-c(iicc,memory)
	
	parametersModel<-list(Order=iicc$q,HXmax=iicc$Hmax,correction_1rOrdre=iicc$correction_1rOrdre,Entropy=iicc$HX,D=iicc$D,Mperfil=iicc$Mperfil,interA=iicc$interA,interB=iicc$interB)
	y<-list(model=iicc$Divergence,parameterModel=parametersModel)	
	return(y)
    
}

	
