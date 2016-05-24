BuildTakensVector <-
function(HRVData, Data, m, tau) {
# -------------------------------------
# Calculates Takens expanded vectors
# -------------------------------------

	if (HRVData$Verbose) {
		cat("** Creating Takens expanded vectors **\n")
		cat("   m: ", m, " Tau: ", tau, "\n", sep="")
	}
	
	N = length(Data)
	jump = tau
	maxjump = (m-1)*jump
	jumpsvect = seq(0,maxjump,jump)
	numjumps = length(jumpsvect)
	numelem = N-maxjump
	DataExp = matrix(nrow=numelem,ncol=numjumps)
	
	for (i in 1:numelem) {
		DataExp[i,1:numjumps] = Data[jumpsvect+i]
	}

	return(DataExp)
}

