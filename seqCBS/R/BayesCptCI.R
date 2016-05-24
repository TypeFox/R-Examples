BayesCptCI <-
function(cases, controls, CBSRes, stepSize="adaptive", adaptMaxMix=80, alpha=0.05, epsilon=10^-4, epsCDF=10^-4, verbose=FALSE) {
	timeCITotal = proc.time()[3]
	timeCIPreProcess = proc.time()[3]
	inStepSize = stepSize
	
	nCas = length(cases)
	nCon = length(controls)
	nTotal = nCas + nCon
	p = nCas/nTotal

	if(verbose)	{
		print(paste("nCas:", nCas, " nCon:", nCon, " nTotal:", nTotal))
	}
	combCasCon = CombineCaseControlC(cases, controls)
	combL = combCasCon$combL
	combX = combCasCon$combX
	combZ = combCasCon$combZ
	combXCumSum = cumsum(as.numeric(combX))
	casTauCumSum = cumsum(as.numeric(combZ))
	conTauCumSum = combXCumSum - casTauCumSum
	nL = length(combL)
	tauHat = c(1, sort(CBSRes$statHat[,3:4]), nL)
	
	timeCIWi = 0
	timeCIMix = 0
	mixStruct = 0
	wkRes = 0
	timeCIPreProcess = proc.time()[3] - timeCIPreProcess
	
	## Compute Wk
	temp = proc.time()[3]
	if(verbose) {
		print(paste("stepSize =", stepSize))
	}
	if(length(tauHat) < 3) {
		## Done, no change point called, use usual Bayesian CI for p
		pBayesCI = matrix(rep(qbeta(p=c(alpha/2, 1-alpha/2), shape1=nCas+1, shape2=nCon+1), nL), ncol=nL)
		pBayesCI = rbind(1:nL, 1:nL, combL, combL, pBayesCI)
		rownames(pBayesCI) = c("LIndex", "RIndex", "LLabel", "RLabel", "CIL", "CIU")
		timeCITotal = proc.time()[3] - timeCITotal
		timeCIWi = 0
		timeCIMix = 0
		timeCINR = 0
		timeCIRes = list(timeCITotal=timeCITotal, timeCIPreProcess=timeCIPreProcess, timeCIWi=timeCIWi, timeCIMix=timeCIMix, timeCINR=timeCINR)
		return(list(CIRes=pBayesCI, timeCIRes=timeCIRes))
	}

	# Not the null model, need to do something
	# 1. Find w_k to the sides of each tauHat, except the first and last
	nTauHatReal = length(tauHat) - 2
	wkRes = vector("list", nTauHatReal+2)
	wkRes[[1]] = matrix(c(1,0,0,casTauCumSum[nL],conTauCumSum[nL],1), nrow=6, ncol=1)
	wkRes[[nTauHatReal+2]] = matrix(c(nL,casTauCumSum[nL],conTauCumSum[nL],0,0,1), nrow=6, ncol=1)
	for(j in 2:(nTauHatReal+1)) {
		wkRes[[j]] = 0
		tauL = as.numeric(tauHat[j-1])
		tauM = as.numeric(tauHat[j])
		tauR = as.numeric(tauHat[j+1])
		tauBoundL = (tauM + tauL)/2
		tauBoundR = (tauM + tauR)/2
		
		sTauL = casTauCumSum[tauM]-casTauCumSum[tauL]
		tTauL = combXCumSum[tauM]-combXCumSum[tauL]
		sTauR = casTauCumSum[tauR]-casTauCumSum[tauM]
		tTauR = combXCumSum[tauR]-combXCumSum[tauM]
		logLTauM = -log(tTauL+1)-log(tTauR+1)-sum(lchoose(c(tTauL, tTauR), c(sTauL, sTauR)))
		continueAdapt = TRUE
		if(inStepSize=="adaptive") {
			useAdapt = TRUE
			stepSize = 5
		}
		else {
			useAdapt = FALSE
		}
		
		while (continueAdapt) {
			wkRes[[j]] = matrix(c(tauM, sTauL, tTauL-sTauL, sTauR, tTauR-sTauR, 1), nrow=6, ncol=1)
			## compute Wi to the left of tauM
			oldK = tauM
			k = tauM - 1
			wkCur = 1
			sTauLCur = sTauL
			tTauLCur = tTauL
			sTauRCur = sTauR
			tTauRCur = tTauR
			while(k > tauBoundL && wkCur > epsilon) {
				sTauChange = casTauCumSum[oldK]-casTauCumSum[k]
				tTauChange = combXCumSum[oldK]-combXCumSum[k]
				sTauLCur = sTauLCur - sTauChange
				tTauLCur = tTauLCur - tTauChange
				sTauRCur = sTauRCur + sTauChange
				tTauRCur = tTauRCur + tTauChange
				logLTauCur = -log(tTauLCur+1)-log(tTauRCur+1)-sum(lchoose(c(tTauLCur, tTauRCur), c(sTauLCur, sTauRCur)))
				wkCur = exp(logLTauCur - logLTauM)
				if(wkCur > epsilon) {
					wkRes[[j]] = cbind(c(k, sTauLCur, tTauLCur-sTauLCur, sTauRCur, tTauRCur-sTauRCur, wkCur), wkRes[[j]])
				}
				oldK = k
				k = k - stepSize
			}
		
			## compute Wi to the right of tauM
			oldK = tauM
			k = tauM + 1
			wkCur = 1
			sTauLCur = sTauL
			tTauLCur = tTauL
			sTauRCur = sTauR
			tTauRCur = tTauR
			while(k < tauBoundR && wkCur > epsilon) {
				sTauChange = casTauCumSum[k]-casTauCumSum[oldK]
				tTauChange = combXCumSum[k]-combXCumSum[oldK]
				sTauLCur = sTauLCur + sTauChange
				tTauLCur = tTauLCur + tTauChange
				sTauRCur = sTauRCur - sTauChange
				tTauRCur = tTauRCur - tTauChange
				logLTauCur = -log(tTauLCur+1)-log(tTauRCur+1)-sum(lchoose(c(tTauLCur, tTauRCur), c(sTauLCur, sTauRCur)))
				wkCur = exp(logLTauCur - logLTauM)
				if(wkCur > epsilon) {
					wkRes[[j]] = cbind(wkRes[[j]], c(k, sTauLCur, tTauLCur-sTauLCur, sTauRCur, tTauRCur-sTauRCur, wkCur))
				}
				oldK = k
				k = k + stepSize
			}
			if(useAdapt && (ncol(wkRes[[j]]) > adaptMaxMix)) {
				continueAdapt = TRUE
				stepSize = stepSize + 5
			}
			else {
				continueAdapt = FALSE
			}
		}
		if(verbose) {
			print(paste("j=", j, "stepSize=", stepSize, "mixLen=", ncol(wkRes[[j]]), "continueAdapt=", continueAdapt))
		}
	}
	timeCIWi = timeCIWi + proc.time()[3] - temp
	wkLen = sapply(wkRes, ncol)
	estMixLen = 100
	for(i in 1:length(wkLen)) {
		estMixLen = c(estMixLen, as.numeric(sapply(1:wkLen[i], function(x) {(wkLen[i]-x)*wkLen[max(1, i-1)] + x*wkLen[min(i+1, length(wkLen))]})))
	}
	if(verbose) {
		print(paste("StepSize:", stepSize))
		print(quantile(estMixLen, probs=c(1, 0.95, 0.9, 0.5)))
	}

	## 2. For each point (or segment if stepSize > 1), we compute
	##    beta parameters and weight for each mixture component
	temp = proc.time()[3]
	nWkEachTauHat = sapply(wkRes, ncol)
	uniquePoint = unlist(sapply(wkRes, function(x) {return(x[1,])}))
	uniquePoint = rbind(uniquePoint[1:(length(uniquePoint)-1)], uniquePoint[2:length(uniquePoint)])
	nUniquePoint = ncol(uniquePoint)
	if(verbose)	print(nUniquePoint)
	
	mixStruct = vector("list", nUniquePoint)
	counter = 1
	for(j in 2:(nTauHatReal+1)) {
		for(k in 1:nWkEachTauHat[j]) {
			casParamL = as.numeric(sapply(casTauCumSum[wkRes[[j]][1,k:nWkEachTauHat[j]]], function(yy) {yy-casTauCumSum[wkRes[[j-1]][1,]]}))
			conParamL = as.numeric(sapply(conTauCumSum[wkRes[[j]][1,k:nWkEachTauHat[j]]], function(yy) {yy-conTauCumSum[wkRes[[j-1]][1,]]}))
			wkL = as.numeric(sapply(wkRes[[j]][6,k:nWkEachTauHat[j]], function(yy) {yy*wkRes[[j-1]][6,]}))
			casParamR = as.numeric(sapply(casTauCumSum[wkRes[[j]][1,0:(k-1)]], function(yy) {casTauCumSum[wkRes[[j+1]][1,]]-yy}))
			conParamR = as.numeric(sapply(conTauCumSum[wkRes[[j]][1,0:(k-1)]], function(yy) {conTauCumSum[wkRes[[j+1]][1,]]-yy}))
			wkR = as.numeric(sapply(wkRes[[j]][6,0:(k-1)], function(yy) {yy*wkRes[[j+1]][6,]}))
			wkLR = c(wkL, wkR)
			wkLR = wkLR/sum(wkLR)
			mixStruct[[counter]] = matrix(rbind(c(casParamL, casParamR), c(conParamL, conParamR), wkLR), ncol=length(wkLR))
			counter = counter+1
		}
	}
	### Now for the very last one
	mixStruct[[counter]] = matrix(wkRes[[nTauHatReal+1]][4:6,], ncol=ncol(wkRes[[nTauHatReal+1]]))
	mixStruct[[counter]][3,] = mixStruct[[counter]][3,]/sum(mixStruct[[counter]][3,])
	timeCIMix = timeCIMix + proc.time()[3] - temp
	if(verbose) {
		print(summary(sapply(mixStruct, ncol)))
		print(length(mixStruct))
		print("mixStruct Computed")
	}

	## 3. Use Newton-Raphson to find quantile of Beta Mixtures
	##    This is implemented in a C routine for speed
	timeCINR = proc.time()[3]
	CIRes = matrix(0, nrow=6, ncol=nUniquePoint)
	CIRes[1,] = uniquePoint[1,]
	CIRes[2,] = uniquePoint[2,]
	CIRes[3,] = combL[CIRes[1,]]
	CIRes[4,] = combL[CIRes[2,]]
	rownames(CIRes) = c("LIndex", "RIndex", "LLabel", "RLabel", "CIL", "CIU")
	for(i in 1:nUniquePoint) {
		CIRes[5:6,i] = .Call("BayesCptCICompC", as.numeric(mixStruct[[i]][1,]), as.numeric(mixStruct[[i]][2,]), as.numeric(mixStruct[[i]][3,]), as.numeric(alpha), as.numeric(epsCDF), PACKAGE="seqCBS")
		if(verbose) {
			reportAt = seq(from=0.1, to=1, by=0.1)
			if(any(i == floor(nUniquePoint*reportAt))) {
				print(paste("CI Completion Progress at: ", 100*reportAt[i == floor(nUniquePoint*reportAt)], "%", sep=""))
			}
		}
	}
	timeCINR = proc.time()[3] - timeCINR
		
	timeCITotal = proc.time()[3] - timeCITotal
	timeCIRes = list(timeCITotal=timeCITotal, timeCIPreProcess=timeCIPreProcess, timeCIWi=timeCIWi, timeCIMix=timeCIMix, timeCINR=timeCINR)
	return(list(CIRes=CIRes, wkRes = wkRes, mixStruct=mixStruct, timeCIRes=timeCIRes))
}		
