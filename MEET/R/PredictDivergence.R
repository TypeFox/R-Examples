PredictDivergence <-function(iicc){

	require("MEET")
	require("seqinr")

    
    write.fasta <- get("write.fasta",pos="package:seqinr")
    read.fasta <- get("read.fasta",pos="package:seqinr")
    
    training.set<-Factortrans<-matriu<-iicc$Transcriptionfactor

	
	iicc$D <- iicc$correction_1rOrdre <- iicc$probparella<-NULL
	iicc$correction_1rOrdre_m1<-iicc$Mperfil<-iicc$q<- NULL
	iicc$HXmax<-iicc$HX<-iicc$q<-NULL
	iicc$Divergence<-iicc$interA<-iicc$interB<-iicc$classentropy<-NULL
	
	
    Prob<-iicc$background
    #iicc$missing.fun=TRUE
    q<-iicc$q<-iicc$model$parameterModel$Order
		
	if (q==1) {iicc$classentropy<-"Shannon"
		}else{
		iicc$classentropy<-"Renyi"  
	}
	
    
    k<-length(iicc$DNA[[1]])
    
    out<-lapply(seq(1,length(iicc$DNA),1), function(r){vector(mode="numeric",length=(k-ncol(Factortrans)+1))})
	out_diff<-lapply(seq(1,length(iicc$DNA),1), function(r){})
	
	iicc$D<-iicc$model$parameterModel$D
	iicc$HXmax<-iicc$model$parameterModel$HXmax
	iicc$correctioc_1rOrdre<-iicc$model$parameterModel$correction_1rOrdre
	iicc$Entropy<-iicc$model$parameterModel$HX
	iicc$Mperfil<-iicc$model$parameterModel$Mperfil
	iicc$interA<-iicc$model$parameterModel$interA
	iicc$interB<-iicc$model$parameterModel$interB
	iicc$Divergence<-iicc$model$model
    names(iicc$Divergence)<-as.character(c(1:length(iicc$Divergence)))
    nameDivergence<-as.character(c(1:length(iicc$Divergence)))
    lengthname<-length(nameDivergence)

	Results<-NULL
   	ncoltraining<-ncol(iicc$Transcriptionfactor)
    MachineDouble<-.Machine$double.eps
	
	for(i in c(1:length(iicc$DNA))){
		
		index<-results<-resultat<-threshold<-a<-Pvalor<-Index<-NULL
		validation.set_x <- iicc$DNA[[i]]
	    iicc$longDNA<-length(iicc$DNA[[i]])

        timein<-proc.time()
        out[[i]]<-.Call("loopDIVERGENCE",
                                        validationsetx=validation.set_x,
                                        lengthDNA=k,
                                        ncoltraining=ncoltraining,
                                        lengthname=lengthname,
                                        interA=iicc$interA,
                                        interB=iicc$interB,
                                        nameDivergence=nameDivergence,
                                        MachineDouble=MachineDouble,
                                        Divergence=iicc$Divergence,
                                        D=iicc$D,
                                        Mperfil=iicc$Mperfil)
        timeout<-proc.time()
        cat("Run time Divergence Detection",timeout[3]-timein[3],"\n")
     
		 a<-rev(sort(as.vector(t(out[[i]]))))
     	 out_diff[[i]]<-as.vector(t(out[[i]]))
    	  jj<-1
      	 threshold<-pvalue(a[jj],out_diff[[i]])
      
     	 while (threshold<iicc$threshold){     
      		 index<-as.numeric(which(out_diff[[i]]==a[jj]))
      		 if (length(index)==1){
      		 Pvalor<-cbind(Pvalor,threshold)	
      		 Index<-cbind(Index,index)
      		 jj<-jj+1
      		 }else{
      		 	for (ii in c(1:length(index))){
      		 	Pvalor<-cbind(Pvalor,threshold)
      		 	Index<-cbind(Index,index[ii])
      		 		}
      		 	jj<-jj+length(index)
      		 	}
      		 threshold<-pvalue(a[jj],out_diff[[i]])      		}
      		 
     
       results<-matrix(,length(Index),3,dimnames=list(c(1:length(Index)),c("Position","Value","Direction")))
	 
    	  results[,1]<-as.numeric(Index[1,])
    	  results[,2]<-as.numeric(round(as.numeric(Pvalor[1,]),7))

      
      if (iicc$direction!="b"){
	  results[,3]<-iicc$direction
	  }else{
	  if (i==1) {results[,3]<-"f"
	  }else{
	  results[,3]<-"r"
	  }
	}
	
      Sequence<-vector(mode="character", length=length(Index))
      for(ii in c(1:length(Index))){
      	if (results[ii,3]=="f"){j<-1
		 	}else{
		 		j<-2}
		 	

	  sequencia<-t(as.matrix(iicc$DNA[[j]][c(Index[ii]:(Index[ii]+ncol(iicc$Transcriptionfactor)))]))
	  Sequence[ii]<-paste(sequencia, sep="", collapse="")
	  }
      resultat <- cbind(results, Sequence)
      Results<-rbind(Results,resultat)
      }
Results

}

