S.STpiPS<-function(S, x, nh)
{
  S<-as.factor(S)
  S<-as.factor(as.integer(S))
  res<-matrix(NA, nrow = sum(nh), ncol=2)
  pk<-matrix(NA,sum(nh))
  cum<-cumsum(nh)
  
  for(k in 1: length(nh)){
    h <- which(S==k)
    res.h <- S.piPS(nh[k], x[h])
    sam.h <- res.h[,1]
    pik.h <- res.h[,2]
    if(k==1){
      res[1:nh[k],1]<-h[sam.h]
      res[1:nh[k],2]<-pik.h
    }
    if(k>1){
      res[(cum[k-1]+1):(cum[k]),1]<-h[sam.h]
      res[(cum[k-1]+1):(cum[k]),2]<-pik.h
    }
  }
  res
}