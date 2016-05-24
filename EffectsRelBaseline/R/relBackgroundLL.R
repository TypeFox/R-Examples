# 
# Compute log-likelihood test statistic of response values being drawn
# from separate populations for each group, different from background
# values.

relBackgroundLL<-function(resp,groups,backMean,backVar) {
  n<-length(resp)
	
  grMeans<-tapply(resp,groups,mean)
  grNums<-tapply(resp,groups,length)
  
  SSC<-sum(grNums*(grMeans-backMean)^2)
	
	groupMeansSpread<-grMeans[as.factor(groups)]
	diffs<-resp-groupMeansSpread
	SSE<-sum(diffs*diffs)
	
	n/2*(1 + log(SSE) - log(n) - log(backVar)) - (SSE+SSC)/(2*backVar)
}
