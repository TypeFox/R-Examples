mouseGaussGenerator<-function(listSizeOfData,EarDistance){
  dimensionOfData<-2
  dataSeat<-matrix(,,dimensionOfData)
  label<-c()
  #glowa
  n<-listSizeOfData[1]
  r<-2;
  dataTemp=matrix(rnorm(n*dimensionOfData,0,r),n,dimensionOfData);  
  m1=0;
  m2=0;
  mean=matrix(c(rep(m1,n),rep(m2,n)),n,dimensionOfData);
  dataSeat<-rbind(dataSeat,dataTemp+mean)
  label<-c(label,rep(1,n))
  #prawe ucho
  r<-1;
  n<-listSizeOfData[2]
  dataTemp=matrix(rnorm(n*dimensionOfData,0,r),n,dimensionOfData);  
  m1=-1-2-EarDistance;
  m2=1+2+EarDistance;
  mean=matrix(c(rep(m1,n),rep(m2,n)),n,dimensionOfData);
  dataSeat<-rbind(dataSeat,dataTemp+mean)
  label<-c(label,rep(2,n))
  #lewe ucho
  r<-1;
  n<-listSizeOfData[3]
  dataTemp=matrix(rnorm(n*dimensionOfData,0,r),n,dimensionOfData);  
  m1=1+2+EarDistance;
  m2=1+2+EarDistance;
  mean=matrix(c(rep(m1,n),rep(m2,n)),n,dimensionOfData);
  dataSeat<-rbind(dataSeat,dataTemp+mean)
  label<-c(label,rep(3,n))
  
  dataSeat<-tail(dataSeat,n=-1)
  list(data=dataSeat,label=label)
}
#listSizeOfData<-c(1000,500,500)
#EarDistance<-2
#test<-mouseGaussGenerator(listSizeOfData,EarDistance)
#plot(test$data,pch=20,col=test$label)

#write.table(test$data,file="C:\\Users\\admin\\Dropbox\\CEC_plugin_R\\TESTY\\mouse_1\\data.txt",row.names=FALSE,col.names=FALSE);
#write.table(test$label,file="C:\\Users\\admin\\Dropbox\\CEC_plugin_R\\TESTY\\mouse_1\\data_cluster.txt",row.names=FALSE,col.names=FALSE,sep = "\n");

