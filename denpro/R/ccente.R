ccente<-function(levels,items,mass){
#Calculates centers from a collection of level sets.
#center is 1st moment didided by volume.
#
#levels is tasolkm*N-matrix of 1:s and 0:s
#items is N*(2*d)-matrix
#mass is tasolkm-vector
#
#returns N*d-matrix of 1st moments.
#
N<-length(levels[,1])
d<-length(items[1,])/2
res<-matrix(0,N,d)
if (dim(t(levels))[1]==1) tasolkm<-1 else tasolkm<-length(levels[,1]) 
for (i in 1:tasolkm){
  lev2<-change(levels[i,])
  m<-length(lev2)
  vol<-matrix(0,d,1)
  for (j in 1:m){
    ind<-lev2[j]
    rec<-items[ind,]
    vol<-vol+cenone(rec)
  }
  res[i,]<-vol/mass[i]
}
return(t(res))
}

 

