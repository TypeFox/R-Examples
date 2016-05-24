kmeansMask<-function(x){
  resamp<-quantile(x,probs = seq(0,1,length.out = 1000))
  mod<-kmeans(resamp,2)
  thresh<-resamp[which.max(abs(diff(mod$cluster)))]
  mask<-ifelse(x<=thresh,0,1)
  return(mask)
}
