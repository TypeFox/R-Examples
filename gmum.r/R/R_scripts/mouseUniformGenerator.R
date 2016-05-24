mouseGaussGenerator<-function(listSizeOfData,EarDistance){
  dataSeat<-matrix(,,dimensionOfData)
  label<-c()
  dimensionOfData<-2
  #glowa
  n<-listSizeOfData[1]
  r<-2;
  X <- runif(4*n,min=-r,max=r)  # generowanie z rozk³adu jednostajnego na [-1,1]
  Y <- runif(4*n,min=-r,max=r)
  Accept <- X^2+Y^2<r^2    # wektor logiczny
  X <- X[Accept]
  Y <- Y[Accept]
  dataTemp=matrix(c(X,Y),length(X),dimensionOfData);
  dataTemp<-tail(dataTemp,n=n)
  
  m1=0;
  m2=0;
  mean=matrix(c(rep(m1,n),rep(m2,n)),n,dimensionOfData);
  dataSeat<-rbind(dataSeat,(dataTemp+mean))  
  label<-c(label,rep(1,n))
  #prawe ucho
  r<-1;
  n<-listSizeOfData[2]

  X <- runif(4*n,min=-r,max=r)  # generowanie z rozk³adu jednostajnego na [-1,1]
  Y <- runif(4*n,min=-r,max=r)
  Accept <- X^2+Y^2<r^2    # wektor logiczny
  X <- X[Accept]
  Y <- Y[Accept]
  dataTemp=matrix(c(X,Y),length(X),dimensionOfData);
  dataTemp<-tail(dataTemp,n=n)
  
  m1=(-1-2-EarDistance)/sqrt(2);
  m2=(1+2+EarDistance)/sqrt(2);
  mean=matrix(c(rep(m1,n),rep(m2,n)),n,dimensionOfData);
  dataSeat<-rbind(dataSeat,dataTemp+mean)
  label<-c(label,rep(2,n))
  #lewe ucho
  r<-1;
  n<-listSizeOfData[3]
  X <- runif(4*n,min=-r,max=r)  # generowanie z rozk³adu jednostajnego na [-1,1]
  Y <- runif(4*n,min=-r,max=r)
  Accept <- X^2+Y^2<r^2    # wektor logiczny
  X <- X[Accept]
  Y <- Y[Accept]
  dataTemp=matrix(c(X,Y),length(X),dimensionOfData);
  dataTemp<-tail(dataTemp,n=n)
  m1=(1+2+EarDistance)/sqrt(2);
  m2=(1+2+EarDistance)/sqrt(2);
  mean=matrix(c(rep(m1,n),rep(m2,n)),n,dimensionOfData);
  dataSeat<-rbind(dataSeat,dataTemp+mean)
  label<-c(label,rep(3,n))
  
  dataSeat<-tail(dataSeat,n=-1)
  list(data=dataSeat,label=label)
}
listSizeOfData<-c(3000,1000,1000)
EarDistance<--0.5
test<-mouseGaussGenerator(listSizeOfData,EarDistance)
plot(test$data,pch=20,col=test$label)

#write.table(test$data,file="C:\\Users\\admin\\Dropbox\\CEC_plugin_R\\TESTY\\mouse_1\\data.txt",row.names=FALSE,col.names=FALSE);
#write.table(test$label,file="C:\\Users\\admin\\Dropbox\\CEC_plugin_R\\TESTY\\mouse_1\\data_cluster.txt",row.names=FALSE,col.names=FALSE,sep = "\n");


