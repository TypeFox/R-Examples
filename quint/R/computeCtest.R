computeCtest <-
function(pmat,dmats3,w){
    #for the test set dmats3 has only 1 row: i.e., the assignment of the training set
    #compute value of criterion
    selp1<-dmats3==1
    selp2<-dmats3==2
    weight<-pmat[,1]+pmat[,2]
    pmat<-cbind(pmat,weight)
    
     #sum of the effect sizes of the regions belonging to p1 weighted by the cardinalities
    dif1<- ifelse(sum(pmat[selp1,4])==0,0,(t(pmat[selp1,4])%*%pmat[selp1,3])/sum(pmat[selp1,4]))
    
    #difference in treatment outcome in leafs assigned to 1
     toleaf1<-pmat[selp1,3]
    
    #sum of the effect sizes of the regions belonging to p2 weighted by the cardinalities
    dif2<- ifelse(sum(pmat[selp2,4])==0,0, (t(pmat[selp2,4])%*%-pmat[selp2,3])/sum(pmat[selp2,4]) )
     #difference in treatment outcome in leafs assigned to 2
         toleaf2<-pmat[selp2,3]
      ##dif1 or dif2 might be negative in testset!!
     #checkdif1<- sapply(1:length(toleaf1),function(kk,toleaf1){ifelse(toleaf1[kk]<0,0,1)  },toleaf1=toleaf1)
     #checkdif2<- sapply(1:length(toleaf2),function(kk, toleaf2){ifelse( toleaf2[kk]>0,0,1)  }, toleaf2= toleaf2)
     checkdif1<- ifelse(dif1<0,0,1)  
     checkdif2<- ifelse(dif2<0,0,1)  
      if(dif1<0){  dif1<-0 }
     if(dif2<0){  dif2<-0  }
      crit1<- w[1]*log(1+dif1)+w[1]*log(1+dif2)
    crit2<- ifelse(sum(pmat[selp1,4])==0|sum(pmat[selp2,4])==0,NA,w[2]*log(sum(pmat[selp1,4]))+w[2]*log(sum(pmat[selp2,4])) )
    crittot<-crit1+crit2
    critch1<-crit1
     if(checkdif1==0){critch1<-1}
     if(checkdif2==0){critch1<-2}
     if(checkdif1==1&checkdif2==1){critch1<-0}
     if(checkdif1==0&checkdif2==1){ 
     #warning("pooled difference in treatment outcome in P1 (i.e., D1) is negative in testset; the difference is set to zero",call.=FALSE,immediate.=FALSE)
     critch1<-1}
     if(checkdif1==1&checkdif2==0){
     #warning("pooled difference in treatment outcome in P2 (i.e., D2) is negative in testset; the difference is set to zero",call.=FALSE,immediate.=FALSE)
     critch1<-2}
     if(checkdif1==0&checkdif2==0){critch1<-3
      #warning("pooled difference in treatment outcome in P1 and in P2 are of opposite sign in testset; the differences are set to zero",call.=FALSE,immediate.=FALSE)
    }
      if(sum(pmat[selp1,4])!=0&sum(pmat[selp2,4])!=0){critch2<-0}
     if(sum(pmat[selp1,4])==0&sum(pmat[selp2,4])!=0){
     critch2<-1
     #warning("cardinality P1 in test set is zero",call.=FALSE) 
     }
     if(sum(pmat[selp1,4])!=0&sum(pmat[selp2,4])==0){
     critch2<-2
     #warning("cardinality P2 in test set is zero",call.=FALSE) 
     }
     if(sum(pmat[selp1,4])==0&sum(pmat[selp2,4])==0){
     critch2<-3
     #warning("both cardinality P1 and cardinality P2 in test set are zero",call.=FALSE) 
     }
    mat<-as.matrix(cbind(crittot,crit1,critch1,crit2,critch2))
    return(mat)
}
