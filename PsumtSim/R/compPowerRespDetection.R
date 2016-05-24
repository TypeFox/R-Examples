# 
# Compute power to detect responses differing from background.
#
compPowerRespDetection<-function(bkgLevel, respLevel, numCats, numCatsWithResp,
																	normDistribution=FALSE, showProgress=FALSE,
																	numTrialsPerCat=15, numBootIters=1000,
																	numRuns=1000, alpha=0.05) {
	
	if (numCatsWithResp>numCats) {
	  stop("number of categories with response must be less then number of categories")
	}
	
	responseLevels<-c(rep(respLevel,numCatsWithResp),
										 rep(bkgLevel,numCats-numCatsWithResp))
	numSig<-0
	
	for (i in 1:numRuns) {
	  if (showProgress) print(i)
	  
		if (normDistribution) sim1<-simNormCatResp(bkgLevel,responseLevels,numTrialsPerCat)
		else sim1<-simCatResp(bkgLevel,responseLevels,numTrialsPerCat)
		
		# Run the bootstrap test and store result for this run
		sigLevel<-testCatEffectBoot(sim1,numBootIters,testFnc=sumSqCat,backMean=bkgLevel)
		
		if (sigLevel<alpha) numSig<-numSig+1
	}

	numSig
}