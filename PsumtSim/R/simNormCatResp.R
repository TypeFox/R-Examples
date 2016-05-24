#
# Build a simulated set of normally distributed responses 
# to categories of stimuli, returning result as a data frame.
#
simNormCatResp<-function(bkgRate,respRates,numRespsPerCat) {
  N<-length(respRates)*numRespsPerCat
	catLabels<-paste('cat',rep(1:length(respRates),numRespsPerCat),sep='')
	bkg<-rnorm(N,bkgRate,bkgRate)
	resp<-rnorm(N,respRates,respRates)
	data.frame(category=catLabels,bkg=bkg,resp=resp)
}