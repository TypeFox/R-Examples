E.STPPS<-function(y,pk,mh,S){
  S<-as.factor(S)
  y<-cbind(1,y)
  y<-as.data.frame(y)
  names(y)[1] <- "N"
  pk<-as.data.frame(pk)
  
  Strata<-array(NA,c(4,length(mh)+1,dim(y)[2]))
  rownames(Strata)=c("Estimation", "Standard Error","CVE","DEFF")
  colnames(Strata)<-c(levels(S),"Population")
  dimnames(Strata)[[3]]<-names(y)
  S<-as.factor(as.integer(S))
  
  for(k in 1: length(mh)){
    e<-which(S==k)
    ye<-y[e,]
    pke<-pk[e,]
    ye<-as.matrix(ye)
    tye<-matrix(1,1,dim(ye)[1])%*%(ye/pke)/mh[k]
    tye2<-t(matrix(tye,dim(y)[2],mh[k]))
    Vtye<-(1/mh[k])*(1/(mh[k]-1))*colSums((ye/pke-tye2)^2)
    CVe<-100*sqrt(Vtye)/tye
    Nh<-(1/mh[k])*sum(1/pke)
    VMAS<-as.vector((Nh^2)*(1-(mh[k]/Nh))*diag(var(ye))/(mh[k]))
    DEFF<-Vtye/VMAS
    Strata[1,,][k,]<-tye
    Strata[2,,][k,]<-sqrt(Vtye)
    Strata[3,,][k,]<-CVe
    Strata[4,,][k,]<-DEFF
  }
  
  m=sum(mh)
  
  for(i in 1:dim(y)[2]){
    Strata[1,,][(length(mh)+1),][i]<-sum(Strata[,,i][1,][1:length(mh)])
    Strata[2,,][(length(mh)+1),][i]<-sqrt(sum(Strata[,,i][2,][1:length(mh)]))
    Strata[3,,][(length(mh)+1),][i]<-100*sqrt(Strata[2,,][(length(mh)+1),][i])/Strata[1,,][(length(mh)+1),][i]
    N=Strata[1,4,1]
    VMAST<-(N^2)*(1-(m/N))*var(y[,i])/(m)
    Strata[4,,][(length(mh)+1),][i]<-(Strata[2,,][(length(mh)+1),][i])/(VMAST)
  }
  return(Strata)
}