wordmatch<-function(file1,file2,n){
  Dir1<-paste(file1,".csv",sep="")
  Dir2<-paste(file2,".csv",sep="")
  t2<-readLines(Dir1)
  t3<-readLines(Dir2)
  t2<-tolower(t2)
  t3<-tolower(t3)
  t2<-strsplit(t2,",")
  t3<-strsplit(t3,",")
  k<-llply(t2,function(x){llply(t3,function(y) which(x %in% y))})
  k1<-melt(k)
  freq=NULL
  k2<-count(k1,c("L1","L2"))
  k3<-subset(k2,freq>=n)
  rm<-t2[k3$L1]
  rm1<-t3[k3$L2]
  k3<-k3[,-3]
  k3["Pair"]<-NA
  k3["Sentence"]<-NA
  k3$Pair<-rm
  k3$Sentence<-rm1
  k3
}



