CalculateRfromCorrelation <-
function(HRVData, Data, m, tau, Cra, Crb) {
# -------------------------------------
# Calculates ra and rb from Correlation
# -------------------------------------
  warning(
    paste("--- Warning: CalculateRfromCorrelation() is deprecated ---",
          "  --- Use CalculateCorrDim() instead ---",
          "  --- See help for more information!! ---",
          sep="\n")
  )
	randC = matrix(nrow=2, ncol=2)

	DataExp = BuildTakensVector(HRVData,Data,m=m,tau=tau)
  #	numelem = nrow(DataExp)

	if (HRVData$Verbose) {
		cat("** Calculating R from Correlation **\n")
	}

	mutualDistance = dist(DataExp,method="maximum")
#	ra = min(mutualDistance)
#	rb = max(mutualDistance)

	numelem = length(mutualDistance)
	rs=quantile(mutualDistance,probs=c(0.005,0.75))
	ra=rs[1]
	rb=rs[2]

	Cmra = length(mutualDistance[mutualDistance<=ra])/numelem
	Cmrb = length(mutualDistance[mutualDistance<=rb])/numelem

	if (HRVData$Verbose) {
		cat("   ra: ", ra, "\n", sep="")
		cat("   rb: ", rb, "\n", sep="")
	}

	if (HRVData$Verbose) {
		cat("   Cmra: ", Cmra*100, "%\n", sep="")
		cat("   Cmrb: ", Cmrb*100, "%\n", sep="")
	}

	randC[1,1] = ra
	randC[1,2] = rb
	randC[2,1] = Cmra
	randC[2,2] = Cmrb
	
	return(randC)
}

