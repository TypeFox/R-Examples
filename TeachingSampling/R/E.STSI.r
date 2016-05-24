E.STSI<-function(S,Nh,nh,y){
  S<-as.factor(S)
  y<-cbind(1,y)
  y<-as.data.frame(y)
  names(y)[1] <- "N"
  
  Strata<-array(NA,c(4,length(nh)+1,dim(y)[2]))
  rownames(Strata)=c("Estimation", "Standard Error","CVE","DEFF")
  colnames(Strata)<-c(levels(S),"Population")
  dimnames(Strata)[[3]]<-names(y)
  S<-as.factor(as.integer(S))
  
  for(k in 1: length(nh)){
    e<-which(S==k)
    ye<-y[e,]
    ye<-as.matrix(ye)
    tye<-matrix(1,1,dim(ye)[1])%*%(ye*(Nh[k]/nh[k]))
    Vtye<-diag((Nh[k]^2)*(1-(nh[k]/Nh[k]))*var(ye)/(nh[k]))
    CVe<-100*sqrt(Vtye)/tye
    VMAS<-diag((Nh[k]^2)*(1-(nh[k]/Nh[k]))*var(ye)/(nh[k]))
    DEFF<-Vtye/VMAS
    Strata[1,,][k,]<-tye
    Strata[2,,][k,]<-sqrt(Vtye)
    Strata[3,,][k,]<-CVe
    Strata[4,,][k,]<-DEFF
  }
  
  N=sum(Nh)
  n=sum(nh)
  
  for(i in 1:dim(y)[2]){
    Strata[1,,][(length(nh)+1),][i]<-sum(Strata[,,i][1,][1:length(nh)])
    Strata[2,,][(length(nh)+1),][i]<-sqrt(sum(Strata[,,i][2,][1:length(nh)]^2))
    Strata[3,,][(length(nh)+1),][i]<-100*(Strata[2,,][(length(nh)+1),][i])/Strata[1,,][(length(nh)+1),][i]
    VMAST<-(N^2)*(1-(n/N))*var(y[,i])/(n)
    Strata[4,,][(length(nh)+1),][i]<-(Strata[2,,][(length(nh)+1),][i])^2/(VMAST)
  }
  return(Strata)
}