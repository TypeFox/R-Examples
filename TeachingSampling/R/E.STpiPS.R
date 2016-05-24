E.STpiPS<-function(y,pik,S){
  S<-as.factor(S)
  y<-cbind(1,y)
  y<-as.data.frame(y)
  names(y)[1] <- "N"
  pik<-as.data.frame(pik)
  nh <- c(table(S))
  
  Strata<-array(NA,c(4,length(nh)+1,dim(y)[2]))
  rownames(Strata)=c("Estimation", "Standard Error","CVE","DEFF")
  colnames(Strata)<-c(levels(S),"Population")
  dimnames(Strata)[[3]]<-names(y)
  S<-as.factor(as.integer(S))
  
  for(k in 1: length(nh)){
    nhe <- nh[k]
    e<-which(S==k)
    ye<-y[e,]
    pike<-pik[e,]
    ye<-as.matrix(ye)
    tye<-matrix(1,1,dim(ye)[1])%*%(ye/pike)
    #-------------------
    ck <- (1-pike)*(nhe/(nhe-1))
    P1 <- as.matrix(colSums(ck*ye/pike))
    P2 <- sum(ck)
    ystar <- t(P1%*%t(pike/P2))
    P3 <- ck/(pike^2)
    #--------------------
    Vtye<-colSums(P3*((ye-ystar)^2))
    CVe<-100*sqrt(Vtye)/tye
    Nhe<-sum(1/pike)
    VMAS<-as.vector((Nhe^2)*(1-(nhe/Nhe))*diag(var(ye))/(nhe))
    DEFF<-Vtye/VMAS
    Strata[1,,][k,]<-tye
    Strata[2,,][k,]<-sqrt(Vtye)
    Strata[3,,][k,]<-CVe
    Strata[4,,][k,]<-DEFF
  }
  
  for(i in 1:dim(y)[2]){
    Strata[1,,][(length(nh)+1),][i]<-sum(Strata[,,i][1,][1:length(nh)])
    Strata[2,,][(length(nh)+1),][i]<-sqrt(sum(Strata[,,i][2,][1:length(nh)]))
    Strata[3,,][(length(nh)+1),][i]<-100*sqrt(Strata[2,,][(length(nh)+1),][i])/Strata[1,,][(length(nh)+1),][i]
    
    N <- sum(1/pik)
    n <- sum(nh)
    VMAST<-(N^2)*(1-(n/N))*var(y[,i])/(n)
    Strata[4,,][(length(nh)+1),][i]<-(Strata[2,,][(length(nh)+1),][i])/(VMAST)
  }
  return(Strata)
}