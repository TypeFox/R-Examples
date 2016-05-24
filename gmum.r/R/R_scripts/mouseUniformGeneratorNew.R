or=function(l){
  bool=FALSE
  for(i in l){
    bool = (bool || i)
  }
  bool
}
l<-c(TRUE,TRUE,FALSE)
or(l)

mouseGaussGenerator<-function(sizeOfData,earDistance,dimensionOfData){
  dataSeat<-matrix(,,dimensionOfData)
  label<-c()
  dimensionOfData<-2
  
  n<-sizeOfData
  rHead<-2;
  rRightEar<-1;
  rLeftEar<-1;
  r<-2*rHead+2*rRightEar+2*rLeftEar
  
  mRightEar1=(-1-2-earDistance)/sqrt(2);
  mRightEar2=(1+2+earDistance)/sqrt(2);
  
  mLeftEar1=(1+2+earDistance)/sqrt(2);
  mLeftEar2=(1+2+earDistance)/sqrt(2);
  
  X <- runif(100*n,min=-r,max=r)  # generowanie z rozk³adu jednostajnego na [-1,1]
  Y <- runif(100*n,min=-r,max=r)
  Accept1 <- (X^2+Y^2<rHead^2)   # wektor logiczny
  Accept2 <- ((X-mRightEar1)^2+(Y-mRightEar2)^2<rRightEar^2)
  Accept3 <- ((X-mLeftEar1)^2+(Y-mLeftEar2)^2<rLeftEar^2)
  Accept<-matrix(c(Accept1,Accept2,Accept3),length(Accept1),3);
  Accept<-apply(Accept,1,or)
  X <- X[Accept]
  Y <- Y[Accept]
  dataTemp=matrix(c(X,Y),length(X),dimensionOfData);
  dataSeat<-tail(dataTemp,n=-1)
  dataSeat<-dataTemp[1:n,]
  
  list(data=dataSeat)
}
sizeOfData<-3000
dimensionOfData = 2
earDistance<--0.2
test<-mouseGaussGenerator(sizeOfData,earDistance,dimensionOfData)
plot(test$data,pch=20)
length(test$data)
#write.table(test$data,file="C:\\Users\\admin\\Dropbox\\CEC_plugin_R\\TESTY\\mouse_1\\data.txt",row.names=FALSE,col.names=FALSE);
#write.table(test$label,file="C:\\Users\\admin\\Dropbox\\CEC_plugin_R\\TESTY\\mouse_1\\data_cluster.txt",row.names=FALSE,col.names=FALSE,sep = "\n");


