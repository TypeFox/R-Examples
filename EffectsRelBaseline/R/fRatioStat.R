#
# F-ratio test statistic for effect of groups.
#
fRatioStat<-function(resp,groups) {
  
  grMeans<-tapply(resp,groups,mean)
  grNums<-tapply(resp,groups,length)
  grandMean<-mean(resp)
  n<-length(resp)
	
	grpFactor<-as.factor(groups)
	MSTr <- sum( grNums * (grMeans - grandMean)^2)
	diffs <- resp - grMeans[grpFactor]
	MSE = sum(diffs^2)
	
  (MSTr/(nlevels(grpFactor)-1))/(MSE/(n-nlevels(grpFactor)))
}