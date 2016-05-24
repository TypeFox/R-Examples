ScanBIC <-
function(combX, combZ, tauHat, lik0, nTotal) {
	## Modified BIC
	nCpts = length(tauHat)-1
	lik1=0
	nCaseWindow = rep(0, nCpts)
	nTotalWindow = rep(0, nCpts)
	for(i in 1:nCpts) {
		nCaseWindow[i] = sum(combZ[tauHat[i]:tauHat[i+1]])
		nTotalWindow[i] = sum(combX[tauHat[i]:tauHat[i+1]])
		pWindow = nCaseWindow[i]/nTotalWindow[i]
		if(0<pWindow && pWindow < 1) {
			lik1 = lik1+nCaseWindow[i]*log(pWindow)+(nTotalWindow[i]-nCaseWindow[i])*log(1-pWindow)
		}
	}
	BICpenalty = -0.5*sum(log(nTotalWindow)) - nCpts*log(length(combX)) + 0.5*log(nTotal)
	mBIC = lik1-lik0+BICpenalty
	return(mBIC)
}

