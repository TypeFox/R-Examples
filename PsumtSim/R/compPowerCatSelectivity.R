# 
# Compute power to detect responses which are category selective using
# the F-ratio statistic.
#
compPowerCatSelectivity<-function(respRates,
																	normDistribution=FALSE, showProgress=FALSE,
																	numTrialsPerCat=15, numBootIters=1000,
																	numRuns=1000, alpha=0.05) {
	
	numSig<-0
	
	for (i in 1:numRuns) {
	  if (showProgress) print(i)
	  
		if (normDistribution) sim1<-simNormCatResp(0.0,respRates,numTrialsPerCat)
		else sim1<-simCatResp(0.0,respRates,numTrialsPerCat)
		
		# Run the bootstrap test and store result for this run
		sigLevel<-testCatEffectBoot(sim1,numBootIters,testFnc=fRatioStat)
		
		if (sigLevel<alpha) numSig<-numSig+1
	}

	numSig
}