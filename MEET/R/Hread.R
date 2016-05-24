Hread<-function(training.set.mes.rand,val.set,iicc){
	
	HXmax   <- iicc$HXmax
    Herror<-iicc$Herror
	Redundancia_corregida<- iicc$Redundancia_corregida
	p<-vector("numeric",length(val.set))
	p[val.set=="A"]<-1;p[val.set=="T"]<-2
	p[val.set=="C"]<-3;p[val.set=="G"]<-4
    p<-as.character(p)
	names(iicc$Entropy)<-as.character(c(1:length(iicc$Entropy)))
    ncoltraining<-ncol(training.set.mes.rand)
    
    HX<-.Call("readH",iicc$Entropy,p,ncoltraining)
    
	out_det <- diffInstructions(training.set.mes.rand, HX, HXmax,Herror, Redundancia_corregida)
	
	out<-(1/(out_det+.Machine$double.eps))
	}
