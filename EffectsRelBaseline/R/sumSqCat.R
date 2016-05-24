# Weighted sum of squares of distance of mean values of each group
# from the background values.
# Hypothesis is fixed by the background mean and groups.
sumSqCat<-function(resp,groups,backMean) {
  
  grMeans<-tapply(resp,groups,mean)-backMean
  grNums<-tapply(resp,groups,length)
  
  sum(grMeans * grMeans * grNums,na.rm=T)
}
  