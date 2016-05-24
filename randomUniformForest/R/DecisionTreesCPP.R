# <OWNER> = Saip CISS
# <ORGANIZATION> = QueensBridge Quantitative
# <YEAR> = 2014

# LICENSE 
# BSD 3-CLAUSE LICENSE

# Copyright (c) 2014, Saip CISS
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# END OF LICENSE 

CheckSameValuesInLabels <- function(Y, n, nodeMinSize)
{
	flagZeros = TRUE
	if (n > nodeMinSize)
	{	if (sum(Y, na.rm = TRUE) !=  n*Y[1]) { flagZeros = FALSE } 	}
	return(flagZeros)
}

CheckSameValuesInAllAttributes <- function(X, n, m, nodeMinSize)
{
	flagZeros = TRUE
	if (n > nodeMinSize)
	{ 	
		i = 1
		Cte <- sum(X[i,,drop = FALSE],na.rm = TRUE)/m	
		while(flagZeros & (i < n))
		{
			if (Cte != X[i,1])
			{  flagZeros = FALSE  }
			else
			{
				if (Cte != X[i+1,1])
				{	flagZeros = FALSE	}
			}
			i = i + 1
		}			
	}	
	return(flagZeros)
}

options.filter <- function(X, Y, nodeMinSize, o.bootstrap = FALSE, o.treeSubsampleRate = FALSE, o.treeOverSampling = 0, o.targetClass = -1, o.OOB = FALSE, o.treeRebalancedSampling = FALSE) 
{
	dimX = dim(X)
	n = nn = dimX[1]
	p = dimX[2]
	follow.idx = 1:n
	flagBalanced = FALSE
	if (!is.numeric(o.treeRebalancedSampling))
	{
		if (o.treeRebalancedSampling)
		{
			classesObject = table(Y)
			minLengthClasses =  which.min(classesObject)
			classes = as.numeric(names(classesObject))
			sampleSize= classesObject[minLengthClasses]
			classesSample = NULL			
			for (i in 1:length(classes))
			{ 	classesSample = c(classesSample, sample(which(Y == classes[i]), sampleSize)) 	}			
			follow.idx = sortCPP((1:n)[classesSample])
			n = length(follow.idx)			
			if (o.OOB) { OOBBalancedIdx = (1:nn)[-follow.idx]; flagBalanced = TRUE }
		}
	}
	else
	{
		if (max(o.treeRebalancedSampling) < 1)
		{	o.treeRebalancedSampling = round(o.treeRebalancedSampling*n,0) }
		if (o.treeSubsampleRate != 1) { stop("rebalanced sampling is not currently compatible with sub-sampling")  }
		classesObject = table(Y)
		minLengthClasses =  which.min(classesObject)
		classes = as.numeric(names(classesObject))
		classesSample = NULL		
		for (i in 1:length(classes))
		{ 	
			classesSampleSize = which(Y == classes[i])
			classesSample = c(classesSample, sample(classesSampleSize, o.treeRebalancedSampling[i], 
			replace = if (classesSampleSize >  o.treeRebalancedSampling[i]) { FALSE } else { TRUE })) 	
		}	
		follow.idx = sortCPP((1:n)[classesSample])
		n = length(follow.idx)
		if (o.OOB) { OOBBalancedIdx =  (1:nn)[-follow.idx]; flagBalanced = TRUE }
	}	
	if (o.treeOverSampling != 0)
	{ 
		if (o.targetClass == (-1)) { stop ("no class to target") }			
		oversample <- overSampling(Y[follow.idx], which.class = o.targetClass, proportion = o.treeOverSampling)	
		follow.idx = sortCPP(c(oversample$C1, oversample$C2))
		n = length(follow.idx)
		
		if (o.treeOverSampling == -1) { o.treeSubsampleRate = 1 }
	}	
	OOBIdx = 0
	if (o.treeSubsampleRate != 1)
	{	
		follow.idx = sortCPP(follow.idx)
		nn = sort(sample(follow.idx, floor(o.treeSubsampleRate*length(follow.idx)), replace = o.bootstrap))
		if (o.OOB)	{  OOBIdx = follow.idx[-nn]	}
		follow.idx = rmNA(follow.idx[nn])		
		if (o.OOB)
		{
			idx = vector()
			for (i in 1:length(follow.idx))
			{	idx = c(idx, which(follow.idx[i] == OOBIdx)) }		
			if (length(idx) > 0) { OOBIdx = OOBIdx[-idx]	}	
		}
		n = length(follow.idx)
	}	
	if ( (o.bootstrap) & (o.treeSubsampleRate == 1) )
	{  
	    bootstrap.idx = sortCPP(sample(n, n, replace = TRUE))				
		follow2.idx = follow.idx[bootstrap.idx]
		
		while(CheckSameValuesInAllAttributes(X[follow2.idx,], n, p, nodeMinSize) | CheckSameValuesInLabels(Y[follow2.idx], n, nodeMinSize) )
		{	bootstrap.idx = sample(n, n, replace = TRUE); follow2.idx = follow.idx[bootstrap.idx]	}
		if (o.OOB)
		{ 	 
			OOBIdx = follow.idx[1:n][-bootstrap.idx]
			rm.idx = find.idx(follow2.idx, OOBIdx, sorting =FALSE)
			if (length(rm.idx) > 0)  { OOBIdx = OOBIdx[-rm.idx] }			
		}		
		follow.idx = follow2.idx
	}
	if (flagBalanced & o.OOB) {  OOBIdx = unique(c(OOBIdx, OOBBalancedIdx))  }	
	return(list(X = X[follow.idx,], Y = Y[follow.idx], n = n, p = p, OOBIdx = OOBIdx, followIdx = follow.idx))
}

uniformDecisionTree <- function(X, Y, nodeMinSize = 1, maxNodes = Inf, treeFeatures = floor(4/3*ncol(X)), 
treeDepth = Inf, 
treeDepthControl = NULL, 
getSplitAt = "random", 
featureSelectionRule = c("entropy", "gini", "random", "L2", "L1"),
regression = FALSE, 
bootstrap = FALSE, 
treeSubsampleRate = 1, 
treeClasswt = NULL, 
treeOverSampling = 0,
targetClass = -1,  
treeRebalancedSampling = FALSE, 
OOB = FALSE, 
treeBagging = FALSE, 
randomFeature = FALSE, 
treeCatVariables = NULL,
outputPerturbation = FALSE, 
unsupervised = FALSE,
moreThreadsProcessing = FALSE, 
moreNodes = FALSE) 
{
	{
		set.seed(sample(1e9,1))
		gainFunction = featureSelectionRule[1]				
		randomSubspace = FALSE
		if (is.character(treeFeatures))
		{ 
			treeFeatures = sample(2:(2*ncol(X)), 1) 
			randomSubspace = TRUE
		}		
		randomNodesize = FALSE
		if (is.character(nodeMinSize))
		{ 
			nodeMinSize = sample(seq(1,min(50,nrow(X)), by = 5),1)
			randomNodesize = TRUE
		}		
		if (is.character(maxNodes)) { maxNodes = sample(3:floor(nrow(X)),1) }				
		if (treeFeatures == 1) { randomFeature = TRUE }				
		if (outputPerturbation) 
		{  
			minValue = if (treeOverSampling == 0) { 0.05 } else { 1/(length(Y)-1) }
			Y <- outputPerturbationSampling(Y, whichClass = targetClass, sampleSize = max(min(treeOverSampling, 1), minValue), 
			regression = regression)
			treeOverSampling = 0 			 
		}
		if (moreThreadsProcessing) { idxList = list();	kk = 1; nodeVector = vector() }
		if (!moreNodes)
		{
			object.filtered <- options.filter(X, Y, nodeMinSize, o.bootstrap = bootstrap, o.treeSubsampleRate = treeSubsampleRate, o.treeOverSampling = treeOverSampling, o.targetClass = targetClass, o.OOB = OOB, o.treeRebalancedSampling = treeRebalancedSampling)
				
			nX = nrow(X)
			init.Y = rep(0, nX)
			init.X = matrix(0, nX, ncol(X)) 
			init.Y = object.filtered$Y
			init.X = object.filtered$X
			n = object.filtered$n
			p = object.filtered$p
			OOBIdx = object.filtered$OOBIdx
			followIdx = object.filtered$followIdx
		}
		else
		{
			init.Y = Y
			init.X = X
			n = nrow(X)
			p = ncol(X)
			OOBIdx = 0
		}
		if (is.character(treeDepth)) { treeDepth = sample(3:round(log(nX)/log(2),0), 1) }
		if (nodeMinSize >= n) { stop("Minimal number of observations (nodesize) cannot be greater than sample size.\n") }
		m = treeFeatures
		if (treeBagging) { m = min(p, treeFeatures); treeFeatures = m }		
		mCoeff = m/p
		iGain = iGain.tmp = rep(0, m)
		minNodes = 3
		thresholdLimit = 0 
		if (!regression) 
		{ 
			if (gainFunction == "random") {	gainFunction = sample(  c("entropy", "gini"), 1) }
			classes = sortCPP(unique(init.Y))
			nClasses = length(classes) 
			YProb = vector(length = nClasses)			
			if (!is.null(treeDepthControl))
			{	
				YProb = tabulate(Y)
				if (gainFunction == "entropy") { L2Limit  <- crossEntropyCPP(YProb/n) } 
				else {  L2Limit  <- giniCPP(YProb/n) }				
				if (is.character(treeDepthControl)) {  limitCoeff = sample(2:128, 1) }
				else  {  limitCoeff = treeDepthControl	}
			}			
			if ( is.null(treeClasswt) ) { newTreeClasswt <- rep(0, nClasses) }
			else  {  newTreeClasswt = treeClasswt }
		}
		else
		{  
			if (gainFunction == "random") {	gainFunction = sample(c("L1", "L2"), 1) }
			classes = init.Y 
			nClasses = 1		   
			if (!is.null(treeDepthControl))
			{	
				if (gainFunction == "L2") 
				{ 
					L2Limit <- L2DistCPP(sum(init.Y, na.rm = TRUE)/n, init.Y)
					if (is.character(treeDepthControl)) {  limitCoeff = sample(4:16, 1) }
					else  {  limitCoeff = treeDepthControl	}
				} 
				else 
				{ 					
					L2Limit <- L1DistCPP(sum(init.Y, na.rm = TRUE)/n, init.Y) 
					if (is.character(treeDepthControl)) {  limitCoeff = sample(2:4, 1) }
					else  {  limitCoeff = treeDepthControl	}
				}				
			}
		}
		nrowTree = if (treeDepth == Inf) { n } else { 2^treeDepth }
		if (maxNodes != Inf) { nrowTree = maxNodes }
		if (nodeMinSize > 1) { nrowTree = floor(nrowTree/nodeMinSize) + 1 } 
		nodes.nb = new.nodes.nb = averageLeafWeight = vector(length = max(10, floor(log(nrowTree))))
		Tree = new.Tree = matrix(data = Inf, ncol = 7, nrow = nrowTree)
		if (n > 100) { 	dataSpaceList = rmVarList = vector("list", n)	}
		else { dataSpaceList = vector("list", n);  rmVarList = vector("list", 10*n) }
		dataSpaceList[[1]] = (1L:n)
		updated.nodes = 0L	
		endCondition = k = idxBestFeature = allNodes = nCte = 1L
		n.dup = maxNumberOfNodes = n
		prev.iGain = vector()
		rmVar = NULL
	}	
	{
		if (moreNodes == FALSE) { rm(object.filtered) }
		allDim = featuresIdx = 1L:p
		flagOptimization_1 = 0L		
		if (treeBagging)  { flagOptimization_1 = 1L }		
		if (randomFeature | outputPerturbation | treeBagging) { lowCorr = FALSE }
		else
		{
		   lowCorr <- if (sample(0:1, 1) == 1) { TRUE } else { FALSE }
		}		
		cteDepth = exp(treeDepth*log(2))
		cteDepth2 = if (is.null(treeDepthControl)) { Inf }
					else
					{	if (regression) { 1/8*cteDepth 	} else { 1/2*cteDepth  }	}
	}
	while (endCondition)
	{		
		set.seed(sample(1e9,1))
		{
			n = length(dataSpaceList[[k]])									
			nodeMinSizeCandidate = max(floor(n/2),1)
			if ( randomNodesize & (k > nodeMinSizeCandidate) ) { nodeMinSize <- sample(nodeMinSizeCandidate,1)	}			
			Y <- init.Y[dataSpaceList[[k]]] 
			X <- init.X[dataSpaceList[[k]], allDim, drop = FALSE]			
			pureNode = FALSE
			flag = 0
			pp = length(rmVarList[[k]])
		}
		allSameValues <- CheckSameValuesInAllAttributes(X, n, p, nodeMinSize)
		if ( allSameValues | ((k >= minNodes) & ( (pp >= p) | (k >= cteDepth) | (k >= maxNodes))) )
		{	
			if ( ((n != 0) & (dataSpaceList[[k]][1] != 0)) | (pp >= p))
			{  
				Leaf <- leafNode(Y, n, nClasses, classes, l.classwt = treeClasswt, l.regression = regression)				
				if (regression)  { Tree[k,] <- c(genericNode(k, -1), Leaf) }
				else { Tree[k,] <- c(genericNode(k, -1), Leaf$Class) }
				nodes.nb[k] = n
				
				if (moreThreadsProcessing & !allSameValues)
				{ 
					if (length(unique(Y)) > 1)
					{
						idxList[[kk]] = dataSpaceList[[k]]
						nodeVector[kk] = k
						kk = kk + 1	
					}
				}				
			}
			else
			{	
				previous.k <- which(Tree == k, arr.ind = TRUE)[1] 
				previousIdx = dataSpaceList[[previous.k]]
				Y <- init.Y[previousIdx]				
				nn = length(previousIdx)		
				Leaf <- leafNode(Y, nn, nClasses, classes, l.classwt = treeClasswt, l.regression = regression)

				if (moreThreadsProcessing)
				{ 
					allSameValues <- CheckSameValuesInAllAttributes(init.X[previousIdx, ], nn, p, nodeMinSize)
					if ( (length(unique(Y)) > 1) & !allSameValues )
					{
						idxList[[kk]] = previousIdx
						nodeVector[kk] = previous.k
						kk = kk + 1	
					}
				}		
				
				if (regression)  { Tree[k,] <- c(genericNode(k, -1), Leaf) }
				else { Tree[k,] <- c(genericNode(k, -1), Leaf$Class) }
				nodes.nb[k] = nn
			}
			updated.nodes = updated.nodes + 2			
			if (!is.null(treeClasswt) & (!regression)) 	
			{  averageLeafWeight[k] <- sum(treeClasswt*Leaf$nb, na.rm = TRUE)/nodes.nb[k]	}
		}
		else
		{
			if ((k >= minNodes) & CheckSameValuesInLabels(Y, n, nodeMinSize) ) 
			{ 
				Tree[k,] <- c( genericNode(k, -1), leafNode(Y, n, nClasses, classes, l.classwt = treeClasswt, 
							which.values = "output", l.regression = regression) )				
				nodes.nb[k] = n 
				pureNode = TRUE
				if (!is.null(treeClasswt) & (!regression)) 	{  	averageLeafWeight[k] = treeClasswt[Y[1]]	}
			}
			else
			{
				if (!regression)
				{
					YProb = tabulate(Y)
					if (is.null(treeClasswt))
					{
						if (gainFunction == "entropy") { YGain <- crossEntropyCPP(YProb/n) } 
						else { YGain <- giniCPP(YProb/n) }
					}
					else
					{
						nClassesLength = length(YProb)
						if (gainFunction == "entropy") 
						{ YGain <- crossEntropyCPP(treeClasswt[1:nClassesLength]*YProb/n) } 
						else {  YGain <- giniCPP(treeClasswt[1:nClassesLength]*YProb/n) }
					}
				}
				else
				{  
					YGain = if (!is.null(treeDepthControl)) { 0 } else { runif(1,0,0.01)  }
					classes = 0 
				}  
				if ( !is.null(rmVarList[[k]]) )
				{	
					previousFeaturesIdx = featuresIdx
					featuresIdx = allDim[-rmVarList[[k]]]
					if (length(featuresIdx) == 0) { featuresIdx = previousFeaturesIdx }
					if (treeBagging) { m = length(featuresIdx) }
					else  { m = max(1, floor(mCoeff*length(featuresIdx))) }
				}
				else	
				{ 	
					featuresIdx = allDim
					m = treeFeatures	
				}
				if ( ((k > cteDepth2) & (m > 1)) | (!is.null(treeDepthControl) & (m > 1)) ) 
				{   
					if (regression) 
					{ 
						if (p <= 32)
						{	candidateVariables = sample(featuresIdx, 16*p, replace = TRUE) }
						else
						{ 	
							if(n <= 32) { candidateVariables = sample(featuresIdx, 4*p, replace = TRUE) } 
							else { candidateVariables = sample(featuresIdx, m, replace = TRUE) }
						}
					}
					else { candidateVariables = sample(featuresIdx, 2*p, replace = TRUE) }	
				}
				else
				{ 				
					if (!flagOptimization_1 & !randomSubspace & lowCorr)
					{
						randomLimit <- if (regression)  { sample(3:sample(3:7,1),1) }  else { 3 }
											
						if ( (k <= randomLimit) & (m >= 2) )
						{   
							U.coeff = runif(1) 
							candidateVariables = sample(featuresIdx, max(2,floor(U.coeff*4*p)), replace = TRUE)	
						}
						else
						{ flagOptimization_1 = 1 }					
					}
					else
					{ candidateVariables <- sample(featuresIdx, m, replace = !treeBagging) }
				}								
				if (randomSubspace)
				{	
					U.coeff <- 2*runif(1)
					candidateVariables <- sample(featuresIdx, max(2, floor(U.coeff*p)), replace = TRUE)	
				}
				np = n*p
				if (length(candidateVariables) > (np)) { candidateVariables = sample(candidateVariables, np) }
				candidateVariables <- sortCPP(candidateVariables)
				catVar = 0 
				if (!is.null(treeCatVariables))
				{
					catVar = NULL 
					for (i in seq_along(treeCatVariables)) 
					{ 
						matchCatVar <- which(candidateVariables == treeCatVariables[i])
						if (length(matchCatVar) > 0) {  catVar = c(catVar, matchCatVar)  }
					}					
					if (!is.null(catVar)) 
					{ 
						XTmp <- X[, candidateVariables, drop = FALSE]
						catVarOld = catVar
						catVar = candidateVariables[catVar]												
						for (j in 1:length(catVar))
						{
							tempXcat = X[,catVar[j]]
							uniquetempXcat = unique(tempXcat)
							if (length(uniquetempXcat) > 1)
							{
								dummy <- sample(uniquetempXcat, 2)
								dummyIdx <- which(tempXcat == dummy[1])
								tempXcat[dummyIdx] = dummy[1]
								tempXcat[-dummyIdx] = dummy[2]
							}
							XTmp[,catVarOld[j]]  = tempXcat
						}
					}
					else
					{ 
						catVar = 0 
						XTmp <- X[, candidateVariables, drop = FALSE]
					}
				}
				else
				{	XTmp <- X[, candidateVariables, drop = FALSE]	}
				if (regression | unsupervised) { variablesLength <- checkUniqueObsCPP(XTmp[, ,drop = FALSE])	} 
				else  
				{ 
					nX = nrow(X)
					sampleIdx = sample(nX, min(nX,16)) 
					variablesLength <- checkUniqueObsCPP(XTmp[sampleIdx, ,drop = FALSE])
				}
				idxVariableLength <- which(variablesLength != 1)
				nIdxVariableLength <- length(idxVariableLength)
				nCurrentFeatures <- length(variablesLength)
				if (nIdxVariableLength > 0) 
				{		
					uniqueObs = variablesLength[-idxVariableLength]   
					lengthUniqueObs <- length(uniqueObs)
					if (lengthUniqueObs > 0) { rmVar <- sortCPP(unique(candidateVariables[uniqueObs])) }
					XTmp <- XTmp[,idxVariableLength, drop = FALSE] 
					candidateVariables = candidateVariables[idxVariableLength]
					if (regression | unsupervised) { thresholds <- runifMatrixCPP(XTmp[, ,drop = FALSE]) }
					else { thresholds <- runifMatrixCPP(XTmp[sampleIdx, ,drop = FALSE])  }
					if (nIdxVariableLength == 1)
					{
						nCurrentFeatures = 1													
						if (!randomFeature)							
						{ 
							if (regression) 
							{ 
								if (gainFunction == "L2") {	iGain <- YGain - L2InformationGainCPP(Y, XTmp, thresholds) }
								else { iGain <- YGain - L1InformationGainCPP(Y, XTmp, thresholds) }
							}
							else  
							{  
								flagIGFunction = 1
								if (gainFunction != "entropy") { flagIGFunction = -1 }								
								iGain <- YGain - entropyInformationGainCPP(Y, XTmp, thresholds, classes, nClasses, newTreeClasswt, flagIGFunction)
							}
						}										
					}					
					else 
					{
						nCurrentFeatures = nIdxVariableLength									
						if (!randomFeature)
						{
							if (regression)  
							{ 	
								if (gainFunction == "L2") {	iGain <- YGain - L2InformationGainCPP(Y, XTmp, thresholds) }
								else { iGain <- YGain - L1InformationGainCPP(Y, XTmp, thresholds) }								
							}
							else 
							{	
								flagIGFunction = 1
								if (gainFunction != "entropy") { flagIGFunction = -1 }								
								iGain <- YGain - entropyInformationGainCPP(Y, XTmp, thresholds, classes, nClasses, newTreeClasswt, flagIGFunction)
							}
						}
					}
				}
				na.gain = which(is.na(iGain) | (iGain == Inf) | (iGain == -Inf))
				na.gain.length = length(na.gain)											
				if (!randomFeature)
				{
					if (na.gain.length == nCurrentFeatures) {  iGain = rep(0, nCurrentFeatures) }
					else {	if (na.gain.length > 0) {  iGain[na.gain] = 0 }  }
				}				
				if (randomFeature & (nIdxVariableLength > 0))
				{
				    idxBestFeature <- if (length(thresholds) > 1) { sample(seq_along(thresholds), 1) } else  { 1 }							
					if (regression)
					{	
						if (gainFunction == "L2") 
						{ 	
							iGain[idxBestFeature] <- YGain - 
							L2InformationGainCPP(Y, XTmp[,idxBestFeature, drop = FALSE], thresholds[idxBestFeature])
						}
						else
						{
							iGain[idxBestFeature] <- YGain - 
							L1InformationGainCPP(Y, XTmp[,idxBestFeature, drop = FALSE], thresholds[idxBestFeature])
						}
					}
					else  
					{  
						flagIGFunction = 1
						if (gainFunction != "entropy") { flagIGFunction = -1 }								
						iGain[idxBestFeature] <- YGain - 
						entropyInformationGainCPP(Y, XTmp[, idxBestFeature, drop = FALSE], 
						thresholds[idxBestFeature], classes, nClasses, newTreeClasswt, flagIGFunction)
					}
				}
				else
				{	
					if (randomFeature) { idxBestFeature <- sample(candidateVariables,1) }
					else { idxBestFeature <- randomWhichMax(iGain) }
				}					
				iGain <- abs(round(iGain,4))
				if (!is.null(treeDepthControl) & (k > minNodes)) { thresholdLimit = L2Limit/(limitCoeff*2*k)  }
				if  (max(iGain) <= thresholdLimit)
				{	
					Leaf <- leafNode(Y, n, nClasses, classes, l.classwt = treeClasswt, l.regression = regression)					
					if (regression)  { Tree[k,] <- c(genericNode(k, -1), Leaf) }
					else { Tree[k,] <- c(genericNode(k, -1), Leaf$Class) }					
					if (!is.null(treeClasswt) & (!regression)) 
					{ averageLeafWeight[k] = sum(treeClasswt*Leaf$nb, na.rm = TRUE)/n }
					nodes.nb[k] = n
					pureNode = TRUE
				}
			}				
			if (pureNode)
			{	updated.nodes = updated.nodes + 2; iGain = iGain.tmp	}
			else
			{
				# patch for some unusual conditions
				if  (max(iGain) == 0)
				{	
					Leaf <- leafNode(Y, n, nClasses, classes, l.classwt = treeClasswt, l.regression = regression)		
					if (regression) { Tree[k,] <- c(genericNode(k, -1), Leaf) }
					else { Tree[k,] <- c(genericNode(k, -1), Leaf$Class) }				
					if (!is.null(treeClasswt) & (!regression)) 
					{ averageLeafWeight[k] = sum(treeClasswt*Leaf$nb, na.rm = TRUE)/n }					
					nodes.nb[k] = n
					pureNode = TRUE
					updated.nodes = updated.nodes + 2
					iGain = iGain.tmp	
				}
				else
				{
					belowK = 2*k - updated.nodes
					aboveK = belowK + 1
					if (aboveK > (maxNumberOfNodes)) 
					{ 
						rmVarList = append(rmVarList, vector("list", n.dup)) 
						nCte = nCte + 1L
						maxNumberOfNodes = nCte*n.dup
					}				
					if (!is.null(rmVarList[[k]]) & !is.null(rmVar))	
					{	rmVarList[[belowK]] = rmVarList[[aboveK]] =  unique(c(rmVarList[[k]], rmVar)) }					
					bestFeature = candidateVariables[idxBestFeature]
					bestThreshold = round(thresholds[idxBestFeature],4)
					XObsOfBestFeature <- X[,bestFeature]				
					subIdx = idxLow = which(XObsOfBestFeature <= bestThreshold)					
					nLowIdx = length(idxLow)  
					nHighIdx = n - nLowIdx
					Tree[k,] <- fullNode(k, belowK, aboveK, bestFeature, bestThreshold)
					nodes.nb[k] = n
					if ( (nLowIdx == n) | (nLowIdx == 0) )
					{ dataSpaceList[[belowK]] = 0L  } 
					else 
					{						
						dataSpaceList[[belowK]] = dataSpaceList[[k]][subIdx]
						flag = 1
						tempXObs = XObsOfBestFeature[idxLow]
						if ( sum(tempXObs, na.rm = TRUE) == (nLowIdx*tempXObs[1]) )  
						{  rmVarList[[belowK]] = unique(c(rmVarList[[belowK]], bestFeature)) }
					}							
					if ( (nHighIdx == 0) | (nHighIdx == n) )  
					{ dataSpaceList[[aboveK]] = 0L; allNodes = aboveK }
					else 
					{	
						dataSpaceList[[aboveK]] = dataSpaceList[[k]][-subIdx]
						tempXObs = XObsOfBestFeature[-idxLow]
						if ( sum(tempXObs, na.rm = TRUE) == (nHighIdx*tempXObs[1])) 
						{  rmVarList[[aboveK]] = unique(c(rmVarList[[aboveK]], bestFeature)) }
						allNodes = aboveK
					}
				}
			}
			prev.iGain[k] = iGain[idxBestFeature] 
			iGain = iGain.tmp
		}		
		rmVarList[[k]] = 0L
		rmVar = NULL 
		k = k + 1L
		if (k > allNodes) {  endCondition = 0	}								
		if ( (k < minNodes) & (endCondition == 0) )
		{ 
			set.seed(sample(1e9,1))
			if (moreThreadsProcessing) { idxList = list();	kk = 1; nodeVector = vector() }
			iGain = iGain.tmp = rep(0, treeFeatures)
			prev.iGain = vector()
			k = endCondition = idxBestFeature = allNodes = 1L			
			if (n.dup > 100) { 	dataSpaceList = rmVarList = vector("list", n.dup)	}
			else { dataSpaceList = vector("list", n.dup);  rmVarList = vector("list", 10*n.dup) }
			dataSpaceList[[1]] = 1L:n.dup
			updated.nodes = flagOptimization_1 = 0				
			Tree = new.Tree = matrix(data = Inf, ncol = 7, nrow = nrowTree)
			nodes.nb = averageLeafWeight = vector(length = max(10, floor(log(nrowTree))))			
			rmVar = NULL
			allDim = 1L:p
			if (treeBagging)  { flagOptimization_1 = 1 }
			thresholdLimit = 0
			if (!regression) 
			{ 
				if (gainFunction == "random") {	gainFunction = sample(  c("entropy", "gini"), 1) }
				classes = sortCPP(unique(init.Y)) 
				nClasses = length(classes) 
				YProb = vector(length = nClasses)				
				if (!is.null(treeDepthControl))
				{	
					YProb = tabulate(Y)
					if (gainFunction == "entropy") { L2Limit  <- crossEntropyCPP(YProb/n) } 
					else {  L2Limit  <- giniCPP(YProb/n) }					
					if (is.character(treeDepthControl)) {  limitCoeff = sample(2:128, 1) }
					else  {  limitCoeff = treeDepthControl	}
				}
			}
			else
			{  
				if (gainFunction == "random") {	gainFunction = sample(  c("L1", "L2"), 1) }
				classes = init.Y 
				nClasses = 1			   
				if (!is.null(treeDepthControl))
				{	
					if (gainFunction == "L2") 
					{ 
						L2Limit <- L2DistCPP(sum(init.Y, na.rm = TRUE)/n, init.Y)
						if (is.character(treeDepthControl)) {  limitCoeff = sample(4:16, 1) }
						else  {  limitCoeff = treeDepthControl	}
					} 
					else 
					{ 					
						L2Limit <- L1DistCPP(sum(init.Y, na.rm = TRUE)/n, init.Y) 
						if (is.character(treeDepthControl)) {  limitCoeff = sample(2:4, 1) }
						else  {  limitCoeff = treeDepthControl	}
					}				
				}
			}	
		}		
		if ( (dim(Tree)[1] - k) < 2 ) 
		{ 
			Tree = rbind(Tree, new.Tree)
			nodes.nb = c(nodes.nb,new.nodes.nb)
		}
	}
	Tree.size = which((Tree == Inf), arr.ind = TRUE)[1] - 1
	Tree = Tree[1:Tree.size,]		
	if (Tree.size == 1) {  Tree = matrix(Tree, nrow = 1) }
	if (is.null(treeClasswt)) 
	{ 
		Tree = cbind(Tree, nodes.nb[1:Tree.size], prev.iGain[1:Tree.size])
		rownames(Tree) = Tree[,1]
		Tree = Tree[,-1]
		if (regression) 
		{ 
			colnames(Tree) = c("left daughter", "right daughter", "split var", "split point", "status", "prediction", "nodes", "L2Dist") 
		}
		else 
		{	
			colnames(Tree) = c("left daughter", "right daughter", "split var", "split point", "status", "prediction", "nodes", "Gain")	
		}
	}
	else
	{
		if (regression) { Tree = cbind(Tree, nodes.nb[1:Tree.size], prev.iGain[1:Tree.size]) }
		else { 	Tree = cbind(Tree, nodes.nb[1:Tree.size], prev.iGain[1:Tree.size], averageLeafWeight[1:Tree.size]) }
		
		rownames(Tree) = Tree[,1]
		Tree = Tree[,-1]
		if (regression) 
		{ colnames(Tree) = c("left daughter", "right daughter", "split var", "split point", "status", "prediction", "nodes", "L2Dist") }
		else 
		{	
			colnames(Tree) = c("left daughter", "right daughter", "split var", "split point", "status", "prediction", "nodes", "Gain", "avgLeafWeight" ) 
		}
	}
	NApredictionIdx = which(is.na(Tree[,"prediction"]))
	if (length(NApredictionIdx) > 0) 
	{ 
		if (regression) Tree[NApredictionIdx, "prediction"] = sum(init.Y)/length(init.Y)
		else Tree[NApredictionIdx, "prediction"] = sample(init.Y, 1)
	}
	if (OOB) 
	{	
		if (moreThreadsProcessing) 
		{ 
			return(list(Tree = Tree, idxList = idxList, nodeVector = nodeVector, OOB.idx = OOBIdx, 
			followIdx = followIdx)) 
		}
		else 
		{ return(list(Tree = Tree, OOB.idx = OOBIdx))	}
	}
	else 
	{	
		if (moreThreadsProcessing) 	
		{ return(list(Tree = Tree, idxList = idxList, nodeVector = nodeVector, followIdx = followIdx)) }
		else 
		{ return(list(Tree = Tree)) }	
	}	
}

genericNode <- function(k, status) c(k, 0, 0, 0, 0, status)
fullNode <- function(k, left.node, right.node, attribute, split.point)  c(k, left.node, right.node, attribute, split.point, 1, 0) 
leafNode <- function(Y, n, nClasses, classes, which.values = "input", l.regression = FALSE, l.classwt = NULL)
{
	if (which.values == "input")
	{
		if (l.regression) 
		{ 
			if (n == 1) { return(Y) }
			else { 	return(sum(Y, na.rm = TRUE)/n) }		
		}
		else
		{ 	
			classNb = rep(0, nClasses)
			classNb[1] = sum(Y == classes[1], na.rm = TRUE)
			if (nClasses == 2) { classNb[2] =  n - classNb[1]}
			else
			{
				for (i in 2:nClasses) 
				{ classNb[i] = sum(Y == classes[i], na.rm = TRUE) }
			}			
			if (!is.null(l.classwt)) {	classNb = l.classwt*classNb }				
			return (list(Class = randomWhichMax(classNb), nb = classNb))
		}			
	}
	else
	{ return(Y[1]) }
}

reduce.trees <- function(object, X, Y, 
newDepth = NULL, 
newMaxNodes = NULL, 
atRandom = FALSE,
regression = TRUE)
{
	n = nrow(X)
	nObject = nrow(object)
	if (!is.null(newDepth)) {  newMaxNodes = 2^newDepth } 
	if (atRandom) 
	{ 
		set.seed(sample(1e9,1))
		newMaxNodes = sample(floor(0.125*nObject):nObject, 1)
	}
		
	if (newMaxNodes < nObject)
	{
		predictedObject <- predictDecisionTree(object, X)
		newObject = object[1:newMaxNodes, ]
				
		nodeVector = which(newObject[,"left daughter"] > newMaxNodes)
		newObject[nodeVector,c(1:4,8)] = 0
		newObject[nodeVector, "status"] = -1
		logNodeVector = round(log(nodeVector),0)
		nLevels = logNodeVector 
		
		for (i in seq_along(nodeVector))
		{	
			values <- rewind.trees(object, nodeVector[i], X = X, nLevels = nLevels[i])$idx	
			if (regression) { newObject[nodeVector[i], "prediction"] = sum(Y[values])/length(values) }
			else {  newObject[nodeVector[i], "prediction"] = as.numeric(names(which.max(table(Y[values])))) }
		}
		
		nodeVector = which(newObject[,"right daughter"] > newMaxNodes)
		if (length(nodeVector) > 0)
		{
			logNodeVector = round(log(nodeVector),0)
			nLevels = logNodeVector 
			appendedNewObject = matrix(0, length(nodeVector), 8)
			colnames(appendedNewObject) = colnames(newObject)
			appendedNewObject[, "status"] = -1
			
			for (i in seq_along(nodeVector))
			{	
				values <- rewind.trees(object, nodeVector[i], X = X, nLevels = nLevels[i])$idx	
				if (regression) { appendedNewObject [i, "prediction"] = sum(Y[values])/length(values) }
				else { appendedNewObject[i, "prediction"] = as.numeric(names(which.max(table(Y[values])))) }
				appendedNewObject[, "nodes"] = length(values)
			}
			
			newObject = rbind( newObject, appendedNewObject)
			rownames(newObject) = 1:nrow(newObject)
		}
		object = newObject
	}
	
	return(object)
}

rewind.trees <- function(Tree, currentNode, X = NULL, nLevels = 3)
{
	rulesMatrix = matrix(NA, nLevels, 3)
	NAIdx = NULL
	for (i in 1:nLevels)
	{
		idx = which(Tree[, "left daughter"] == currentNode)
		lengthIdx = length(idx)
		if (lengthIdx == 0) 
		{	idx = which(Tree[, "right daughter"] == currentNode); lengthIdx = length(idx) }
		splitVar = Tree[idx, "split var"][1]
		splitPoint = Tree[idx, "split point"][1]
		if (lengthIdx > 0) 	{	rulesMatrix[i,] = c(splitVar, splitPoint, -1)	}
		else  {  NAIdx = c(NAIdx, i) }
		currentNode = idx[1]
	}
	colnames(rulesMatrix) = c("split var", "split point", "direction")
	
	if (!is.null(NAIdx)) 
	{ 
		rulesMatrix = rulesMatrix[-NAIdx, , drop = FALSE] 
		nLevels = nLevels - length(NAIdx)
	}
	
	if (is.null(X))
	{	return(rulesMatrix)	}
	else
	{
		n = nrow(X) 
		XIdx = tempIdx = rep(0, n)
		
		if (nLevels > 0 )
		{
			for (i in 1:nLevels)
			{
				#if (!is.na(rulesMatrix[i,1]))
				{
					if (rulesMatrix[i,3] == 1) 	{  tempIdx = which(X[ ,rulesMatrix[i,1]] > rulesMatrix[i,2]) }
					else {  tempIdx = which(X[ ,rulesMatrix[i,1]] <= rulesMatrix[i,2]) }
					XIdx[tempIdx] = XIdx[tempIdx] + 1
				}
			}
		}
		else
		{  idx = 1:n }
		
		return(list(rules = rulesMatrix, idx = which(XIdx >= nLevels)))
	}
}


onlineClassify <- function(treeObject, X, c.regression = 0)
{
	p = length(X)
	YonlineLearning = X[p]; X = X[-p]
	classifyObject <- classifyMatrixCPP(treeObject, matrix(X, 1, p))
	nodes.idx = 1
	if (c.regression != 0)  
	{ 	
		if ( (YonlineLearning - classifyObject[1])^2  > c.regression )  {  nodes.idx = 0  } 	
	}
	else 
	{ 	
		if (YonlineLearning == classifyObject[1]) 	{	nodes.idx = 0 	} 	
	}
	return(nodes.idx)
}

predictDecisionTree <- function(treeObject, X, p.onlineLearning = FALSE, Y.p.onlineLearning = 0, p.regression = 0) 
{
	n = nrow(X)
	if (!is.matrix(treeObject)) { treeObject = treeObject$Tree	}
	if (p.onlineLearning) 
	{  
		X = cbind(X, Y.p.onlineLearning)
		predictedObject = vector(length = n) 
		for (i in 1:n) 
		{  predictedObject[i] <- onlineClassify(treeObject, X[i,], c.regression = p.regression)	}
	}
	else 
	{  
		predictedObject = matrix(data = 0, ncol = 3, nrow = n)
		colnames(predictedObject) = c("value", "nodes.data", "depth") 
		predictedObject <- classifyMatrixCPP(treeObject, X)	
	}	
	return(predictedObject)
}
# END OF FILE