MIread<-function(training.set,val.set,iicc){
	
	p<-vector("numeric",length(val.set))
    p[val.set=="A"]<-1;p[val.set=="T"]<-2
	p[val.set=="C"]<-3;p[val.set=="G"]<-4
    names(iicc$Divergence)<-as.character(c(1:length(iicc$Divergence)))
	name<-as.character(c(1:length(iicc$Divergence)))
    lengthname<-length(name)
    ncoltraining<-ncol(training.set)
    interA<-iicc$interA
    interB<-iicc$interB

    Dd<-matrix(.Call("MI",iicc$Divergence,name,lengthname,ncoltraining,p,interA,interB),ncoltraining,ncoltraining)

	#Out trace
	preout<-abs(iicc$D-Dd)*iicc$Mperfil
	#diag(preout)<-0
	out<- 1/sum(preout+.Machine$double.eps)

	}	
			
