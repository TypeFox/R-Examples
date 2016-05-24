.index.G1d<-function(d,cl)
{
  d<-as.matrix(d)
  wgss<-0
  bgss<-0
  for(i in 1:(nrow(d)-1))
  for(j in (i+1):nrow(d)){
    if(cl[i]==cl[j]){
      wgss<-wgss+d[i,j]^2
    }
    else{
      bgss<-bgss+d[i,j]^2
    }
  }
  (bgss/(max(cl)-1))/((sum(d^2)/2-bgss)/(length(cl)-max(cl)))
}


#library(clusterSim)
#x<-cluster.Gen(model=7)
#for(i in 2:10){
#print(paste("klasy",i))
#cl<-pam(x$data,i)$clustering
#print(index.G1(as.matrix(x$data),cl))
#print(index.G1d(dist(x$data),cl))
#print("--------------------------")
#}



