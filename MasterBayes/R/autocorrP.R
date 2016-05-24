"autocorrP"<-function(postP){

  nitt<-dim(postP)[2]/2

  Parcomb<-matrix(paste(postP[,seq(1,nitt*2,2)], postP[,seq(2,nitt*2,2)]), dim(postP)[1], dim(postP)[2]/2)

  lag1<-apply(Parcomb, 1, function(x){sum(x[1:(nitt-1)]!=x[1:(nitt-1)+1])})/(nitt-1)
  lag2<-apply(Parcomb, 1, function(x){sum(x[1:(nitt-2)]!=x[1:(nitt-2)+2])})/(nitt-2)
  lag3<-apply(Parcomb, 1, function(x){sum(x[1:(nitt-5)]!=x[1:(nitt-5)+5])})/(nitt-5)
  lag4<-apply(Parcomb, 1, function(x){sum(x[1:(nitt-10)]!=x[1:(nitt-10)+10])})/(nitt-10)
  lag5<-apply(Parcomb, 1, function(x){sum(x[1:(nitt-50)]!=x[1:(nitt-50)+50])})/(nitt-50)
  lag6<-apply(Parcomb, 1, function(x){sum(x[1:(nitt-100)]!=x[1:(nitt-100)+100])})/(nitt-100)
  
  ACP<-as.matrix(c(mean(lag1-lag2), mean(lag1-lag3),mean(lag1-lag4),mean(lag1-lag5), mean(lag1-lag6)))

  rownames(ACP)<-c("Lag 2", "Lag 5", "Lag 10", "Lag 50", "Lag 100")
  colnames(ACP)<-"P"

  ACP
}

