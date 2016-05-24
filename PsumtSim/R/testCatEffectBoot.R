#
# Test for an effect of category using bootstrapping,
# returning the achieved significance level.
#

testCatEffectBoot<-function(sim,R,testFnc,...) {
		
  testVal<-testFnc(sim$resp,sim$category,...)
  bootRes<-boot(sim,catEffectBootAdaptor,R,testFnc=testFnc,...)
  sum(bootRes$t>testVal)/R
}