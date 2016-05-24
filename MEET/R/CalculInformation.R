CalculInformation <-function(matriu, Prob){
  B<-1
	
       number<-matrix(,4,ncol(matriu))
       R<-vector(mode='logical', length=ncol(matriu))
       for (ii in c(1:(ncol(matriu)))) {
               number[1,ii]<-sum(as.numeric(matriu[,ii]=="A"));number[2,ii]<-sum(as.numeric(matriu[,ii]=="T"))
               number[3,ii]<-sum(as.numeric(matriu[,ii]=="C"));number[4,ii]<-sum(as.numeric(matriu[,ii]=="G"))
	       R[ii]<- sum(as.numeric(matriu[,ii]=="-"))
	       
       }
       transnum<-t(number)
       
       logodds<-score<-frequencies<-matrix(,ncol(matriu),4)
              for (i in c(1:4)){
		frequencies[,i] <-(transnum[,i]+(B*Prob[i]))/(nrow(matriu)+B)
               score[,i]<-frequencies[,i]/Prob[i]
	      
	      

       }
  
  logodds<-log2(score)

  informacio <- matrix(,1,nrow(score))
    for(i in c(1:nrow(score))){
      informacio[1,i]<- frequencies[i,]%*%logodds[i,]
    
    }
    
    matriu_pesos <- matrix(,nrow(score),4)
    
    for(i in c(1:nrow(score))){
     matriu_pesos[i,]<- score[i,]*informacio[i]
    
    }
    
    matriu_inf<-cbind(t(informacio), matriu_pesos)
    matriu_inf

}

