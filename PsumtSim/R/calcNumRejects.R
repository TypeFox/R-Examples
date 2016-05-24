# 
# Calculate number of cases rejected for repeated simulation
# of Poisson background and responses grouped into categories.
#
calcNumRejects<-function(bkg,resps,numRespsPerCat,numSims,calcPValFnc,sigLevel=0.05,...) {
  numRejected<-0
  for (i in 1:numSims) {
    sim<-simCatResp(bkg,resps,numRespsPerCat)
    pVal<-calcPValFnc(sim,...)
    if (pVal<sigLevel) numRejected = numRejected + 1
  }
  numRejected
}