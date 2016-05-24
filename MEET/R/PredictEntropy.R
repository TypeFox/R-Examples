PredictEntropy <-function(iicc){
  require("MEET")

  Prob<-iicc$background
  training.set<-Factortrans<-matriu<-iicc$Transcriptionfactor
    #iicc$missing.fun=TRUE
  
  k<-length(iicc$DNA[[1]])
 
  out<-lapply(seq(1,length(iicc$DNA), 1), function(x){
			return(matrix(,1,k-ncol(Factortrans)+1))})
  out_diff<-lapply(seq(1,length(iicc$DNA), 1), function(x){vector(mode="numeric",length=(k-ncol(Factortrans)+1))})

  q<-iicc$q<-iicc$correction_1rOrdre<-iicc$Herror<-iicc$HXmax<-iicc$Redundancia_corregida<-iicc$ErrorHX<-iicc$classentropy<-iicc$Entropy<-NULL
	
  q<-iicc$q<-iicc$model$parameterModel$Order
  
  if (q==1) {iicc$classentropy<-"Shannon"
		}else{
			iicc$classentropy<-"Renyi"  
	}

  #new line
  timein<-proc.time()
  iicc$Herror<-slot(correction.entropy(q,nrow(training.set),1,iicc),"Herror")[nrow(iicc$Transcriptionfactor)]
  timeout<-proc.time()
  cat("run tine correction entropy",timeout[3]-timein[3],"\n")

  iicc$HXmax<-iicc$model$parameterModel$HXmax
  iicc$ErrorHX<- iicc$model$parameterModel$ErrorHX
  iicc$Redundancia_corregida<-iicc$model$parameterModel$Redundancy
    
  iicc$Entropy<-iicc$model$model
  names(iicc$Entropy)<-as.character(c(1:length(iicc$Entropy)))
  nameEntropy<-c("1","2","3","4")
  Results<-NULL
  
  for(i in c(1:length(iicc$DNA))){
      
      index<-results<-resultat<-threshold<-a<-Pvalor<-Index<-NULL
      validation.set_x <-iicc$DNA[[i]]
      
      #load("Time.RData")
      timein<-proc.time()
      
      out[[i]]<-.Call("loopENTROPY",
                        validationsetx=validation.set_x,
                        lengthDNA=k,
                        ncoltraining=ncol(iicc$Transcriptionfactor),
                        Redundanciacorregida=iicc$Redundancia_corregida,
                        iicc=iicc,
                        Entropy=iicc$Entropy,
                        nameEntropy=nameEntropy,
                        MachineDouble=.Machine$double.eps)
      timeout<-proc.time()
      #time<-timeout[3]-timein[3]
      #TimeEntropy<-c(time,TimeEntropy)
      #save("TimeEntropy","TimeDivergence",file="Time.RData")
      cat("Run time Entropy Detection",timeout[3]-timein[3],"\n")
      
     # out[[i]][1,]<-(sapply(X=c(1:(k-ncol(matriu)+1)),
     #	    FUN = function(X, training.set, validation.set_x,iicc) {
	 #	    seq.rand <-validation.set_x[X:(X+ncol(training.set)-1)]
	 #		Hread(training.set=training.set,val.set=seq.rand,iicc)
	 #	    }, training.set=training.set, iicc=iicc, validation.set_x=validation.set_x))
        
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
      		 threshold<-pvalue(a[jj],out_diff[[i]])      			}
      	
      	
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

	  sequencia<-t(as.matrix(iicc$DNA[[1]][c(Index[ii]:(Index[ii]+ncol(iicc$Transcriptionfactor)))]))
	  Sequence[ii]<-paste(sequencia, sep="", collapse="")
	  }
	  
      resultat <- cbind(results, Sequence)
      Results<-rbind(Results,resultat)
     
     }

Results
}

