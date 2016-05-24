cvolum<-function(levels,items){
#Calculates volumes of set of level sets
#
#levels is tasolkm*N-matrix of 1:s and 0:s
#items is N*(2*d)-matrix
#
#returns N-vector of volumes
#
N<-length(levels[,1])
res<-matrix(0,N,1)
if (dim(t(levels))[1]==1) tasolkm<-1 else tasolkm<-length(levels[,1]) 
for (i in 1:tasolkm){
  lev2<-change(levels[i,])
  m<-length(lev2)
  vol<-0
  for (j in 1:m){
    ind<-lev2[j]
    rec<-items[ind,]
    vol<-vol+massone(rec)
  }
  res[i]<-vol
}
return(t(res))
}


 

