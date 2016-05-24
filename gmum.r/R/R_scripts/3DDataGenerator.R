gaussGenerator<-function(listSizeOfData,dimensionOfData){
  dataSeat<-matrix(,,dimensionOfData)
  label<-c()
  for(i in (1:length(listSizeOfData))){
    n<-listSizeOfData[i]; 
    rotaion<-matrix(rnorm(dimensionOfData*dimensionOfData,0,1), nrow = dimensionOfData) 
    dataTemp=matrix(rnorm(n*dimensionOfData,0,1),n,dimensionOfData);  
    m1=2*rnorm(1,0,1);
    m2=m1;
    mean=matrix(c(rep(m1,n),rep(m2,n)),n,dimensionOfData);
    dataSeat<-rbind(dataSeat,dataTemp%*%rotaion+mean)
    label<-c(label,rep(i,n))
  }
  dataSeat<-tail(dataSeat,n=-1)
  list(data=dataSeat,label=label)
}

listSizeOfData<-c(700,300,1000)
dimensionOfData<-3

test<-gaussGenerator(listSizeOfData,dimensionOfData)

plot(test$data,pch=20,col=test$label)

#write.table(test$data,file="C:\\Users\\admin\\Dropbox\\CEC_plugin_R\\TESTY\\mouse_1\\data.txt",row.names=FALSE,col.names=FALSE);
#write.table(test$label,file="C:\\Users\\admin\\Dropbox\\CEC_plugin_R\\TESTY\\mouse_1\\data_cluster.txt",row.names=FALSE,col.names=FALSE,sep = "\n");

