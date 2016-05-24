CalculPSSM <-function(matriu,Prob){
      B<-1
       number<-matrix(,4,ncol(matriu))
	  
       for (ii in c(1:(ncol(matriu)))) {
               number[1,ii]<-sum(as.numeric(matriu[,ii]=="A"));number[2,ii]<-sum(as.numeric(matriu[,ii]=="T"))
               number[3,ii]<-sum(as.numeric(matriu[,ii]=="C"));number[4,ii]<-sum(as.numeric(matriu[,ii]=="G"))
	       
       }
      
       transnum<-t(number)
       
       logodds<-score<-matrix(,ncol(matriu),4)
       
       for (i in c(1:4)){
               score[,i]<-(((transnum[,i]+(B*Prob[i]))/(nrow(matriu)+B)))/Prob[i]

       }
  
logodds<-log2(score)
logodds
}

