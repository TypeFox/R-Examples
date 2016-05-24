#
# log-likelihood ratio for normal 1-way anova
# A test of category responses being different from the average response.
#
LL1WayAnova<-function(resp,groups) {
  n<-length(resp)

  meanResp<-mean(resp)
  grMeans<-tapply(resp,groups,mean)
  grNums<-tapply(resp,groups,length)
  
  SSRC<-sum(grNums * (grMeans-meanResp)^2)
  
  groupMeansSpread<-grMeans[as.factor(groups)]
	diffs<-resp-groupMeansSpread
	SSE<-sum(diffs*diffs)

  n/2*(log(SSE)-log(SSE+SSRC))
}