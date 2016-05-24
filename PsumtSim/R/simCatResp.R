#
# Build a simulated set of responses to categories of stimuli, 
# returning result as a data frame.
#
simCatResp<-function(bkgRate,respRates,numRespsPerCat) {
  N<-length(respRates)*numRespsPerCat
	catLabels<-paste('cat',rep(1:length(respRates),numRespsPerCat),sep='')
	bkg<-rpois(N,bkgRate)
	resp<-rpois(N,respRates)
	data.frame(category=catLabels,bkg=bkg,resp=resp)
}