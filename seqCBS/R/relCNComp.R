relCNComp <-
function(combX, combZ, tauHatInd, p, alpha) {
	baselineRelCN = p/(1-p)
	nTau = length(tauHatInd)-1
	relCN = rep(0, nTau)
	relGainLoss = matrix(0, nrow=nTau, ncol=2)
	colnames(relGainLoss) = c("testInd", "pval")
	for(i in 1:nTau) {
		tolRead = sum(combX[tauHatInd[i]:tauHatInd[i+1]])
		caseRead = sum(combZ[tauHatInd[i]:tauHatInd[i+1]])
		relCN[i] = (caseRead/(tolRead-caseRead+1))/baselineRelCN
		gainP = (binom.test(caseRead, tolRead, p=p, alternative="greater")$p.value)*(nTau*2)
		lossP = (binom.test(caseRead, tolRead, p=p, alternative="less")$p.value)*(nTau*2)
		if(gainP < alpha) {
			relGainLoss[i,] = c(1, gainP)
		} else if(lossP < alpha) {
			relGainLoss[i,] = c(-1, lossP)
		} else {
			relGainLoss[i,] = c(0, (binom.test(caseRead, tolRead, p=p)$p.value)*(nTau*2))
		}
	}
	return(list(relCN=relCN, relGainLoss=relGainLoss))
}

