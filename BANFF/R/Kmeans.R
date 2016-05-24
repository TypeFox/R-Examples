#####Using Keams for Providing Initial Values
Kmeans=function(rstat)
{
  wholeindex=kmeans(rstat,2)$cluster
  mean1=mean(rstat[which(wholeindex==1)])
  mean2=mean(rstat[which(wholeindex==2)])
  if ((mean2)>(mean1)){
    wholeindex<-wholeindex-rep(1,length(rstat))
  }else{
    wholeindex[which(wholeindex==2)]<-0}
  return(wholeindex)
}
