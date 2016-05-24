# apply continous haar fitz transform
chft <-function(y, dyadicL, intermidL)
{
	res = cht(y, dyadicL, intermidL)
	scalingCoeff  = res$scalingCoeff 
	waveletCoeff = res$waveletCoeff 
	scalingCoeff  = scalingCoeff  + as.numeric(scalingCoeff ==0)
	res = waveletCoeff/sqrt(scalingCoeff )
	return(res)
}

# y: observed input
# dyadicL: number of dydaic scales
# intermidL: number of intermidiate scales
# returns objetcs scalingCoeff and waveletCoeff
# waveletCoeff continuous wavelet transform coefficients
# scalingCoeff   scaling coefficients at each scale
cht <-function(y, dyadicL, intermidL)
{
    waveletTrasform = matrix(0, nrow = dyadicL, ncol = length(y))
	waveletTrasform[1,] = y 
	for (level in 2:dyadicL)
	{
		waveletTrasform[level, ] = waveletTrasform[(level-1),] + shift.sequence(waveletTrasform[(level-1),], -1*(2^(level-2)))
	}
	
	jump = 0
	L = 0
	while (L < intermidL)
	{
		jump = jump+1
		L = 2^(jump+1)-jump-2	
	}
	
	jump = min(jump,dyadicL-2)
	intermidL = 2^(jump+1)-jump-2
	totalL = dyadicL+intermidL
	inbetween = c(rep(1, (dyadicL-jump)),2^(1:jump))
	dyadicscales = cumsum(inbetween)

	waveletCoeff = matrix(0, nrow = max(dyadicscales), ncol = ncol(waveletTrasform))
	scalingCoeff = waveletCoeff
	for (level in 1:dyadicL)
	{
		x = shift.sequence(waveletTrasform[level,],-1*(2^(level-1)))
		waveletCoeff[dyadicscales[level],] = waveletTrasform[level,] - x
		scalingCoeff[dyadicscales[level],] = waveletTrasform[level,] + x
	}

	for (level in (dyadicL-jump):(dyadicL-1))
	{
		waveletAtLevel = waveletTrasform[level,]
		sWaveletAtLevel = shift.sequence(waveletTrasform[(dyadicL-jump-1),], -1*(2^(level-1)))
		for (iLevel in 1:(2^(level-dyadicL+jump+1)-1))
		{
			waveletAtLevel = waveletAtLevel+sWaveletAtLevel;
			x = shift.sequence(waveletAtLevel,-1*(2^(level-1)+iLevel*2^(dyadicL-jump-2)))
			waveletCoeff[(dyadicscales[level]+iLevel),] = waveletAtLevel - x
			scalingCoeff[(dyadicscales[level]+iLevel),] = waveletAtLevel + x
			sWaveletAtLevel = shift.sequence(sWaveletAtLevel,-1*(2^(dyadicL-jump-2)))
		}	
	}
	return(list("waveletCoeff" = waveletCoeff, "scalingCoeff" = scalingCoeff))
}

# shifts wavelet coefficients (for centralizing the transform)
shiftWaveletCoeff<-function(wavletCoeff, dyadicL, intermidL)
{
	totalL = nrow(wavletCoeff)
	jump = 0
	L = 0
	while (L < intermidL)
	{
	   jump = jump+1
	   L = 2^(jump+1)-jump-2
	}
	jump = min(jump,dyadicL-2)
	intermidL = 2^(jump+1)-jump-2
	currShift = 0.5
	for (level in 2:(dyadicL-jump)) # dydaic shifts
	{
		currShift = currShift*2
		wavletCoeff[level,] = shift.sequence(wavletCoeff[level,],currShift);
	}
	d = currShift
	for (level in (dyadicL-jump+1):totalL)
	{
		currShift = currShift+d
		wavletCoeff[level,] = shift.sequence(wavletCoeff[level,],currShift) 
	}
	return(wavletCoeff)
}

# performs convulation (R convolve has accuracy problems)
myConv<-function(u,v)
{
	m = length(u)
	n = length(v)
	size = m+n-1
	w = vector(length = size, mode = "numeric")
	for (k in 1:size)
	{
		jMin = max(1, (k+1-n))
		jMax = min(k,m)
		w[k] = sum(u[jMin:jMax]*v[(k-jMin+1):(k-jMax+1)])
	}
	return(w)
}

convolveIter <-function(f0, iter)
{
	f1 = f0
	for (i in 2:iter)
	{
	   f2 = rep(0, (length(f0)*2-1))
	   f2[seq(1, length(f2), 2)] = f0 # every second place 
	   f1 = myConv(f2,f1)
	   f0 = f2
	}
	return(f1)
}


findWaveletMaximaByScale <- function(wavletCoeff)
{
	
	# first make a a filter to eliminate local maxima that arise due to noise
	f0 = c(.000892313668,    -.001629492013,   -.007346166328, 
	.016068943964,   .026682300156,   -.081266699680, 
	-.056077313316,  .415308407030,   .782238930920,
	.434386056491,   -.066627474263,  -.096220442034, 
	.039334427123,   .025082261845,   -.015211731527,
	-.005658286686,  .003751436157,   .001266561929,
	-.000589020757,  -.000259974552,  .000062339034,
	.000031229876,   -.000003259680,  -.000001784985)
	f0Rev = rev(f0)
	iter = 6
	
	localFilter = myConv(convolveIter(f0Rev,iter), convolveIter(f0, iter))
	localFilter = localFilter/(2^iter)
	
	
	M = nrow(wavletCoeff)
	N = ncol(wavletCoeff)
	maxima = matrix(0, nrow = M, ncol = N)
	for (level in 1:M)
	{
	   waveletCoeffAtScale = wavletCoeff[level,]
	   currMax = max(localFilter)
	   maxPos = which(currMax == localFilter)[1]
	   smootheWaveletCoef = myConv(localFilter, waveletCoeffAtScale)
	   smootheWaveletCoef = smootheWaveletCoef[maxPos:(maxPos+N-1)]
	   # find local maxima 
	   f = abs(smootheWaveletCoef)	   
	   n = length(f);
	   x = min(f)
	   ff = c((x-1), f, (x-1))
	   sLocalMaxima = which(ff[2:(n+1)] >= ff[1:n] & ( ff[2:(n+1)] >= ff[3:(n+2)]))

	   # match maxima in smoothed to precised location
	   pLocalMaxima = c(1, sLocalMaxima, N) 
	   pLocalMaxima = floor((pLocalMaxima[1:(length(pLocalMaxima)-1)] + pLocalMaxima[2:length(pLocalMaxima)])/2)
	   for (i in 1:length(sLocalMaxima))
	   {
		  i1 = pLocalMaxima[i]
		  i2 = pLocalMaxima[i+1]
		  x = wavletCoeff[level,i1:i2]*sign(wavletCoeff[level,sLocalMaxima[i]])
		  sLocalMaxima[i] = which(x ==  max(x))[1]
		  sLocalMaxima[i] = sLocalMaxima[i] - 1 + pLocalMaxima[i]
	   }
	   maxima[level,sLocalMaxima] = 256;
	}
	return(maxima)
}

expandMatrix <-function(m, addRow, addCol)
{
	M = nrow(m)
	N = ncol(m)
	e = matrix(0, nrow = (M + addRow) ,ncol = (N+addCol))
	e[1:M, 1:N] = m
	return(e)
}

getOverlap <-function(lineMinPos,lineMaxPos,lineMinScale,lineMaxScale)
{
	# encourage overlap in terms of location rather than scales
	tempLineMinPos = 2*lineMinPos - lineMaxPos+1
	lineMaxPos = 2*lineMaxPos - lineMinPos-1
	lineMinPos = tempLineMinPos;

	n = length(lineMinPos);
	overlap = matrix(0, nrow = n, ncol = n)
	for(i in 2:n)
	{
	   noOverlaps = (lineMinPos[i] > lineMaxPos[1:i-1]) | (lineMaxPos[i] < lineMinPos[1:i-1])
	   overlapInScale = min(c(lineMaxScale[i]-lineMinScale[1:i-1], lineMaxScale[1:i-1]-lineMinScale[i]))
	   noLocalOverlapsInScale = lineMaxScale[i] - lineMinScale[i] - overlapInScale
	   noGlobalLocalOverlapsInScale = lineMaxScale[1:i-1] - lineMinScale[1:i-1] - overlapInScale
	   dominantOverlap = (overlapInScale > (noLocalOverlapsInScale/4)) | (overlapInScale > (noGlobalLocalOverlapsInScale/4))
	   m = which(!(noOverlaps) & (!(dominantOverlap))) 
	   overlap[i,m] = 1
	   overlap[m,i] = 1
	}
	return(overlap)
}


# preprocessMaxima: merge close maxima, removes unreasonable connections and then merges again
# the output are candidate change points
preprocessMaxima <-function(maxima, wavletCoeff, dyadicL ,intermidL, threshold)
{

	# 1. merge close maxima that occur due to noise 
	# num of observations (left & right), for each scale, 
	# involved in the computation of a balanced Haar wavelet coefficients 
	numObservarionsByScale = 1:(dyadicL + intermidL)
	totalNumOfObs = length(numObservarionsByScale); 
	seqSize = ncol(wavletCoeff);
	for (level in 1:totalNumOfObs)
	{
	   n = numObservarionsByScale[level]
	   for (i in 1:seqSize)
	   {
			if (maxima[level,i] != 0)
			{
				leftObsNum = max(1,(i-n)) 
				righObsNum = min(seqSize,(i+n))
				if (abs(wavletCoeff[level,i]) < max(abs(wavletCoeff[level,leftObsNum:righObsNum])))
				{
					maxima[level,i] = 0
				}
			}
	   
	   }
	}
	
	#2. connect maxima lines across succesive scales
	numScales = nrow(maxima)
	maximaAtCoarseLevel = which(maxima[numScales,] != 0)
	numMaximaAtCoarseLevel = length(maximaAtCoarseLevel)
	maximaPos = matrix(0, nrow = numScales, ncol = numMaximaAtCoarseLevel)
	maximaPos[numScales,] = maximaAtCoarseLevel
	coarsestLines = 1:numMaximaAtCoarseLevel
	
	for (level in (numScales-1):1)
	{
	   maximaAtLevel = which(maxima[level,] != 0)
	   n = length(maximaAtLevel)
	  
	   if (n > 0)
	   {
		   currMaximaLines = rep(1, n) 
		   closestLineInCoarse = vector(length = n, mode = "numeric")
		   
		   for (i in 1:n)
		   {
				x = abs(maximaAtLevel[i]-maximaAtCoarseLevel)
				closestLineInCoarse[i] = which(min(x) == x)[1]
		   }
		   
		   mergedMaximaLines = vector(mode = "numeric")
		   for (j in 1:length(maximaAtCoarseLevel))
		   {
				x = abs(maximaAtCoarseLevel[j]-maximaAtLevel)
				closestPos = which(min(x) == x)
				
				if (closestLineInCoarse[closestPos] == j) # merge
				{
					 currMaximaLines[closestPos] = 0
					 lineNum = coarsestLines[j]
					 if (lineNum > ncol(maximaPos))
					 {
						maximaPos = expandMatrix(maximaPos, 0, (lineNum - ncol(maximaPos)))
					 }
					 mergedMaximaLines[closestPos] = lineNum
					 maximaPos[level,lineNum] = maximaAtLevel[closestPos]	 
				}
			}
			
			mergedMaximaLines[which(is.na(mergedMaximaLines))] = 0
			
			k = ncol(maximaPos)
			indices = which(currMaximaLines != 0)
			currMaximaPos = maximaPos
			for (closestPos in indices)
			{
			  k = k+1
			  maximaPos = expandMatrix(maximaPos, 0, (k - ncol(maximaPos)))
			  mergedMaximaLines[closestPos] = k
			  maximaPos[level,k] = maximaAtLevel[closestPos]
			}
		    coarsestLines = mergedMaximaLines
		    maximaAtCoarseLevel = maximaAtLevel
			   
		}
		  
	}
	
	#3. remove local lines that are not significant 
	lastColWithChange = 0;
	maximaLowLines = matrix(0, nrow = nrow(maximaPos), ncol = ncol(maximaPos))
	for (i in 1:ncol(maximaPos)) 
	{
	   levelsWithMax = which(maximaPos[,i] != 0) # all levels were the maxima of col i is present
	   n = length(levelsWithMax)
	   if (n > 0)
	   {
		   indices = maximaPos[levelsWithMax,i]
		   waveletCoeffAtIndex = vector(length = length(levelsWithMax), mode = "numeric")
		   for (j in 1:length(levelsWithMax))
		   {
			  waveletCoeffAtIndex[j] = wavletCoeff[levelsWithMax[j],indices[j]]
		   }
		   if (max(abs(waveletCoeffAtIndex)) > threshold)
		   {
			  # keep lines of local maxima if at least one maximum has a coeff above/equal to the threshold
			  lastColWithChange = lastColWithChange+1;
			  maximaLowLines[,lastColWithChange] = maximaPos[,i]
		   }
		}
	}
	maximaPos = maximaLowLines[,1:lastColWithChange]
	
	#4. remove unreasonable connections
	nLines = ncol(maximaPos)
	nScales  = nrow(maximaPos) 
	n = nLines
	if (length(n) == 0) # no change points
	{
		maximaPos = NULL
	}else{
		for (i in 1:n)
		{
		   levelsWithMax = which(maximaPos[,i] != 0)
		   indices = maximaPos[levelsWithMax,i]
		   r = length(indices)
		   jumps = 0
		   if (r > 1)
		   {
			  jumps = (indices[2:r] - indices[1:(r-1)])/(levelsWithMax[2:r] - levelsWithMax[1:(r-1)])
		   }
		   threshJumps = sqrt(sum((jumps-mean(jumps))^2)/length(jumps))
		   threshJumps = min(20,floor(3*threshJumps))
		   
		   while(threshJumps %in% abs(jumps))
		   {
				threshJumps = threshJumps+1
		   }
		   meanjumps = mean(jumps)
		   if (length(jumps) == 1)
		   {	   
				meanjumps = 0
				threshJumps = 20
		   }
		   outliers = which(abs(jumps-meanjumps) > threshJumps)
		   if (length(outliers) > 0 )
		   {
			   outliers = sort(outliers)
			   for (j in 1:length(outliers))
			   {
				  outlier = outliers[j]
				  l = levelsWithMax[outlier]
				  nLines = nLines+1
				  maximaPos = expandMatrix(maximaPos,0,1)
				  maximaPos[1:l,nLines] = maximaPos[1:l,i]
				  maximaPos[1:l,i] = 0
			   }
		   }
	 
		}
		
		#5. merge the lines 
		nLines = ncol(maximaPos)
		nScales  = nrow(maximaPos) 
		pmax = max(max(maximaPos)) - min(min(maximaPos[maximaPos > 0]))+1
		lineDensity = rowSums(maximaPos != 0)/pmax
		realLines = rep(1,nLines)
		lineMinPos = vector(length = nLines, mode = "numeric") 
		lineMaxPos = vector(length = nLines, mode = "numeric") 
		lineMinScale = vector(length = nLines, mode = "numeric") 
		lineMaxScale = vector(length = nLines, mode = "numeric")
		nLinePoints = vector(length = nLines, mode = "numeric")	
		#find the minimum and maximum position/scale of each maxima line
		for (i in 1:nLines)
		{
			levelsWithMax = which(maximaPos[,i] != 0)
			currMaximaVals = maximaPos[levelsWithMax,i]
			lineMinPos[i] = min(currMaximaVals)
			lineMaxPos[i] = max(currMaximaVals)
			lineMinScale[i] = min(levelsWithMax)
			lineMaxScale[i] = max(levelsWithMax)
			nLinePoints[i] = length(levelsWithMax)
		}
		lineMedianPos = (lineMaxPos + lineMinPos)/2
		
		# test for ovrelap between values
		overlap = getOverlap(lineMinPos,lineMaxPos,lineMinScale,lineMaxScale)
		#link
		lineLength = lineMaxScale - lineMinScale + 1
		res = sort(-nLinePoints*lineMaxScale, index.return = T)
		sortedLength = res$x * -1
		sortedLines = res$ix
		n = max(which(sortedLength > 1))
		for (i in 1:n)
		{
		   currLine = sortedLines[i]
		   candidates = which(overlap[currLine,] != 0)
		   while (length(candidates) != 0)
		   {
			  jump = Inf
			  lineLinks = matrix(0, nrow = length(candidates),ncol = 4)		 
			 # find end position of bridges 
			 for (j in 1:length(candidates))
			 {		 
				 candidate = candidates[j]
				 if (lineMinScale[candidate] < lineMinScale[currLine])
				 {
					lineLinks[j,2] = lineMinScale[currLine]
					lineLinks[j,4] = min(lineMaxScale[candidate],lineMinScale[currLine])
				 }else{
					lineLinks[j,2] = lineMaxScale[currLine]
					lineLinks[j,4] = max(lineMinScale[candidate],lineMaxScale[currLine])
				 }
				 lineLinks[j,1] = maximaPos[lineLinks[j,2],currLine]
				 lineLinks[j,3] = maximaPos[lineLinks[j,4],candidate]
			  }
			  
			  # test if we can link lines
			  for (j in 1:length(candidates))
			  {
				 p1 = min(lineLinks[j,1],lineLinks[j,3])
				 p2 = max(lineLinks[j,1],lineLinks[j,3])
				 k1 = min(lineLinks[j,2],lineLinks[j,4])
				 k2 = max(lineLinks[j,2],lineLinks[j,4])
				 for (k in k1:k2)
				 {
					if (sum((maximaPos[k,] < p2) & (maximaPos[k,] > p1)) > 0)
					{
						lineLinks[j,1] = Inf
					}
				 }
				 
			  }

			  # when the link across scales is too long, don't use it
			  gap1 = abs(lineLinks[,4] - lineLinks[,2])
			  gap2 = abs(lineLinks[,3] - lineLinks[,1])
			  currLength = lineLength[currLine]
			  linkLengths = gap1/currLength + gap2*lineDensity[lineLinks[,2]]
			  linkIndex = which(gap1 > (lineLength[currLine]/2))
			  if (length(linkIndex > 0))
			  {
				linkLengths[linkIndex] = Inf
			  }
			  # in case two bridges have the same length, let the median position of
			  # the involved lines of local maxima decide
			  linkLengths = linkLengths + abs(lineMedianPos[currLine] - lineMedianPos[candidates])*lineDensity[lineLinks[,2]]
			  minLength = min(linkLengths)
			  linkIndex = which(linkLengths == minLength)[1]
			  if (minLength < Inf)
			  {
				 candidate = candidates[linkIndex]
				 minScale = lineMinScale[currLine]
				 maxScale = lineMaxScale[currLine]
				 prevMaxLine = maximaPos[minScale:maxScale, currLine]
				 maximaPos[1:minScale,currLine] = maximaPos[1:minScale,candidate]
				 maximaPos[maxScale:nScales,currLine] = maximaPos[maxScale:nScales,candidate]
				 maximaPos[minScale:maxScale,currLine] = prevMaxLine							
				 levelsWithMax = which(maximaPos[,currLine]!=0)
				 if (length(levelsWithMax) > 0)
				 {
					 currMaximaVals = maximaPos[levelsWithMax,currLine]
					 lineMinPos[currLine] = min(currMaximaVals)
					 lineMaxPos[currLine] = max(currMaximaVals)
					 lineMinScale[currLine] = min(levelsWithMax)
					 lineMaxScale[currLine] = max(levelsWithMax)
					 maximaPos[,candidate] = 0
					 lineMinPos[candidate] = 0
					 lineMaxPos[candidate] = 0
					 lineMinScale[candidate] = 0
					 lineMaxScale[candidate] = 0
					 realLines[candidate] = 0
					 
					 allCandidates = c(candidates, which(overlap[candidate,] != 0))
					 overlap[allCandidates,allCandidates] = getOverlap(lineMinPos[allCandidates],lineMaxPos[allCandidates],lineMinScale[allCandidates],lineMaxScale[allCandidates])
					 overlap[,candidate] = 0
					 overlap[candidate,] = 0
					 candidates = allCandidates[which(overlap[currLine,allCandidates] != 0)]
				 }
			  }else{
				candidates = c()
			  }
		   }
		}
		x = which(realLines!=0)
		if (length(x) > 0)
		{
			maximaPos = maximaPos[,x]
		}
	}
	
	return(maximaPos)
}



findScaleConfidenceRange<-function(y, position, rangeStart, rangeEnd, scale)
{

	N = length(y); # y are the original obersvations
	leftScale = min(scale,position);
	rightScale = min(scale,N-position);
	maxCoeff = 0
	maxLeftScale = leftScale
	maxRightScale = rightScale
	minLeftPos = 1
	maxLeftPos = position
	minRighttPos = 1
	maxRightPos = N-position
	stepLeft = ceiling((maxLeftPos-minLeftPos)/512)
	stepRight = ceiling((maxRightPos-minRighttPos)/512)
	

	# calculate left and right boundaries
	while (stepLeft > 0 | stepRight > 0)
	{
	   if (stepLeft > 0 & maxLeftPos >= minLeftPos)
	   {
		   mySteps = seq(minLeftPos, maxLeftPos, stepLeft)
		   for (leftScale in mySteps)
		   {
			   localSum = sum(y[(position-leftScale+1):(position+rightScale)])
			   if (localSum > 0 & (rightScale*leftScale >0))
			   {
				 coeff = (mean(y[(position+1):(position+rightScale)])-mean(y[(position-leftScale+1):position]))/sqrt(localSum)*sqrt(leftScale*rightScale)
				 if (abs(coeff) > abs(maxCoeff))
				 {
					maxLeftScale = leftScale
					maxCoeff = coeff
				 }
			  }
		   }
	   }
	   leftScale = maxLeftScale
	   if (stepRight > 0 & maxRightPos >= minRighttPos)
	   {
		   mySteps = seq(minRighttPos, maxRightPos, stepRight)
		   for (rightScale in mySteps)
		   {
			  localSum = sum(y[(position-leftScale+1):(position+rightScale)])
			  if ((localSum > 0) & (rightScale*leftScale > 0))
			  {
					coeff = (mean(y[(position+1):(position+rightScale)])-mean(y[(position-leftScale+1):position]))/sqrt(localSum)*sqrt(leftScale*rightScale)
					if (abs(coeff) > abs(maxCoeff))
					{
						maxRightScale = rightScale
						maxCoeff = coeff
					}
			  }
		   }
	   }
	   rightScale = maxRightScale
	   if (stepLeft > 0 & maxLeftPos >= minLeftPos)
	   {
		   mySteps = seq(minLeftPos, maxLeftPos, stepLeft)
		   for (leftScale in mySteps)
		   {
			  localSum = sum(y[(position-leftScale+1):(position+rightScale)])
			  if ((localSum > 0) & (rightScale*leftScale > 0))
			  {
				coeff = (mean(y[(position+1):(position+rightScale)])-mean(y[(position-leftScale+1):position]))/sqrt(localSum)*sqrt(leftScale*rightScale)
				 if (abs(coeff) > abs(maxCoeff))
				 {
					maxLeftScale = leftScale
					maxCoeff = coeff
				 }
			  
			  }
			}
		}
		
	   leftScale = maxLeftScale
	   minLeftPos = max(leftScale-stepLeft,minLeftPos)
	   maxLeftPos = min(leftScale+stepLeft,maxLeftPos)
	   minRighttPos = max(rightScale-stepRight,minRighttPos)
	   maxRightPos = min(rightScale+stepRight,maxRightPos)
	   stepLeft = max(0,min(stepLeft-1,ceiling((maxLeftPos-minLeftPos)/512)))
	   stepRight = max(0,min(stepRight-1,ceiling((maxRightPos-minRighttPos)/512)))
	}

	# adjust location according to boundaries 
	leftBound = position-leftScale+1
	rightBound = position+rightScale
	rangeStart = max(rangeStart,leftBound)
	rangeEnd = min(rangeEnd,rightBound)
	steps = ceiling((rangeEnd-rangeStart)/512)
	while (steps > 0 & (rangeEnd >= rangeStart) )
	{
	   currMaxPos = position
	   leftBound = position-leftScale+1
	   rightBound = position+rightScale
	   localSum = sum(y[leftBound:rightBound])
	   mySteps = seq(rangeStart, rangeEnd, steps)
	   for (position in mySteps)
	   {
		  leftScale = position-leftBound+1
		  rightScale = rightBound-position
		  if (leftScale > 0 & rightScale > 0)
		  {
			 coeff = (mean(y[(position+1):rightBound])-mean(y[leftBound:position]))/sqrt(localSum)*sqrt(leftScale*rightScale)
		  }else{
			 coeff = 0
		  }
		  if (abs(coeff) > abs(maxCoeff))
		  {
			 currMaxPos = position
			 maxLeftScale = leftScale
			 maxRightScale = rightScale
			 maxCoeff = coeff
		 }
	   }
	   position = currMaxPos
	   leftScale = maxLeftScale
	   rightScale = maxRightScale

	   rangeStart = max(position-steps,rangeStart)
	   rangeEnd = min(position+steps,rangeEnd)
	   steps = min(steps-1,ceiling((rangeEnd-rangeStart)/512))

	}
	return(list("position" = position, "leftScale" = leftScale, "rightScale" = rightScale, "coeff" = maxCoeff))
}


findOptimalPositionWithScaleConstriant<-function(y, position, constraintLeftScale, constraintRightScale)
{
	
	leftScale = constraintLeftScale
	rightScale = constraintRightScale
	maxLeftScale = constraintLeftScale
	maxRightScale = constraintRightScale
	maxLeftPos = maxLeftScale
	maxRightPos = maxRightScale
	minLeftPos = 1
	minRighttPos = 1
	stepLeft = ceiling((maxLeftPos-minLeftPos)/512)
	stepRight = ceiling((maxRightPos-minRighttPos)/512)
	maxCoeff = 0

	# calculate left and right boundaries
	while (stepLeft > 0 | stepRight > 0)
	{
	   if (stepLeft > 0 & maxLeftPos >= minLeftPos)
	   {
		   mySteps = seq(minLeftPos, maxLeftPos, stepLeft)
		   for (leftScale in mySteps)
		   {
			   localSum = sum(y[(position-leftScale+1):(position+rightScale)])
			   if (localSum > 0 & (rightScale*leftScale >0))
			   {
				 coeff = (mean(y[(position+1):(position+rightScale)])-mean(y[(position-leftScale+1):position]))/sqrt(localSum)*sqrt(leftScale*rightScale)
				 if (abs(coeff) > abs(maxCoeff))
				 {
					maxLeftScale = leftScale
					maxCoeff = coeff
				 }
			  }
		   }
	   }
	   leftScale = maxLeftScale
	   if (stepRight > 0 & maxRightPos >= minRighttPos)
	   {
		   mySteps = seq(minRighttPos, maxRightPos, stepRight)
		   for (rightScale in mySteps)
		   {
			  localSum = sum(y[(position-leftScale+1):(position+rightScale)])
			  if ((localSum > 0) & (rightScale*leftScale > 0))
			  {
					coeff = (mean(y[(position+1):(position+rightScale)])-mean(y[(position-leftScale+1):position]))/sqrt(localSum)*sqrt(leftScale*rightScale)
					if (abs(coeff) > abs(maxCoeff))
					{
						maxRightScale = rightScale
						maxCoeff = coeff
					}
			  }
		   }
	   }
	   rightScale = maxRightScale
	   if (stepLeft > 0 & maxLeftPos >= minLeftPos)
	   {
		   mySteps = seq(minLeftPos, maxLeftPos, stepLeft)
		   for (leftScale in mySteps)
		   {
			  localSum = sum(y[(position-leftScale+1):(position+rightScale)])
			  if ((localSum > 0) & (rightScale*leftScale > 0))
			  {
				coeff = (mean(y[(position+1):(position+rightScale)])-mean(y[(position-leftScale+1):position]))/sqrt(localSum)*sqrt(leftScale*rightScale)
				 if (abs(coeff) > abs(maxCoeff))
				 {
					maxLeftScale = leftScale
					maxCoeff = coeff
				 }
			  
			  }
			}
		}
		
	   leftScale = maxLeftScale
	   minLeftPos = max(leftScale-stepLeft,minLeftPos)
	   maxLeftPos = min(leftScale+stepLeft,maxLeftPos)
	   minRighttPos = max(rightScale-stepRight,minRighttPos)
	   maxRightPos = min(rightScale+stepRight,maxRightPos)
	   stepLeft = max(0,min(stepLeft-1,ceiling((maxLeftPos-minLeftPos)/512)))
	   stepRight = max(0,min(stepRight-1,ceiling((maxRightPos-minRighttPos)/512)))
	}
	return(list("leftScale" = leftScale, "rightScale" = rightScale, "coeff" = maxCoeff))
}

findSigPosScale <- function(wavletCoeff, dyadicL, maximaLines, y)
{
	intermidL = nrow(wavletCoeff) - dyadicL
	numObservarionsByScale = 1:(dyadicL + intermidL)
	nLines = ncol(maximaLines)
	sigRes = matrix(0, nrow = nLines, ncol = 5)
	colnames(sigRes) =  c("scale", "position", "leftScale", "rightScale", "coeff")
	for (maximaLine in 1:nLines)
	{
	   levelsWithMax = which(maximaLines[,maximaLine] != 0)
	   lineMaximaLevels = maximaLines[levelsWithMax,maximaLine]
	   lineMaximaCoeff = vector(length = length(levelsWithMax), mode = "numeric")
	   for (i in 1:length(levelsWithMax))
	   {	   
			lineMaximaCoeff[i] = wavletCoeff[levelsWithMax[i],lineMaximaLevels[i]]
       } 
	   i = which(abs(lineMaximaCoeff) ==  max(abs(lineMaximaCoeff)))[1]
	   j = levelsWithMax[i]
	   scale = numObservarionsByScale[j]
	   position = lineMaximaLevels[i]
	   rangeStart = min(lineMaximaLevels)
	   rangeEnd = max(lineMaximaLevels)
	   res = findScaleConfidenceRange(y,position,rangeStart,rangeEnd,scale)
	   sigRes[maximaLine, 1] = j
	   sigRes[maximaLine, 2] = res$position
	   sigRes[maximaLine, 3] = res$leftScale
	   sigRes[maximaLine, 4] = res$rightScale
	   sigRes[maximaLine, 5] = res$coeff
	}
	return(as.data.frame(sigRes))
}
# calculates Akaiki's Information Criterion for a given segmentation (change points)
calcAIC<-function(y, changePoints)
{
	numObs = length(y)
	changePoints = sort(changePoints)
	if (changePoints[1] > 1)
	{
		changePoints = c(1, changePoints)
	}
	numCP = length(changePoints)
	if (changePoints[numCP] < (numObs+1))
	{
		changePoints = c(changePoints, numObs+1)
		numCP = numCP + 1
	}
	expected = vector(length = numObs, mode = "numeric")
	for (p in 1:(numCP-1))
	{
	   p0 = changePoints[p]
	   p1 = changePoints[p+1]
	   expected[p0:(p1-1)] = mean(y[p0:p1-1])
	}

    x = (expected*(as.numeric(expected > 0)) + as.numeric(expected==0))
    logLikelihood = sum(-expected + y*log(x)*(as.numeric(expected>0)) - lfactorial(y))
    indices = which(expected == 0)
    if (length(indices) > 0)
    {
	   if (max(y[indices]) > 0)
	   {
		 logLikelihood = Inf
	   }
    }
	# k = (2*numCP - 3) is number of estimated parameters
	aic = -2*logLikelihood + 2*(2*numCP-3)
	return(aic)					  
}

findOptimizedCP <-function(y, maximaLines, sigVals, eps)
{
	candidateCP1 = c(1, (length(y)+1))
	candidateCP2 = candidateCP1
	aic = calcAIC(y,candidateCP1)
	zVals = rep(1, ncol(maximaLines))
	zValue = max(abs(sigVals$coeff))
	linePos = which(abs(sigVals$coeff) == zValue)[1]
	zVals[linePos] = zValue
	selectedLines = linePos
	sigVals$coeff[linePos] = 0
	remainingLines = which(abs(sigVals$coeff) > eps)
	selectedChangePoints = c()
	while(length(remainingLines) > 0)
	{
	   p = sigVals$position[linePos]
	   candidateCP1 = sort(c(candidateCP1, p))
	   candidateCP2 = c(candidateCP2, p)
	   aic = c(aic, calcAIC(y, candidateCP1))
	   for(currLine in remainingLines)
	   {
		  pos = sigVals$position[currLine]
		  if (pos %in% candidateCP1)
		  {
			 sigVals$coeff[currLine] = 0
		  }else{
			 leftScale = sigVals$leftScale[currLine]
			 rightScale = sigVals$rightScale[currLine]
			 left = pos-leftScale+1
			 right = pos+rightScale
			 if ((pos <  p) & (right > p))
			 {
				rightScale = p-pos
			 }
			 if ((pos > p) & (left < p))
			 {
				leftScale = pos-p+1
			 }
			 
			 res = findOptimalPositionWithScaleConstriant(y,pos,leftScale,rightScale)
			 sigVals$leftScale[currLine] = res$leftScale
			 sigVals$rightScale[currLine] = res$rightScale
			 sigVals$coeff[currLine] = res$coeff
		  }
	   }
	   zValue = max(abs(sigVals$coeff))
	   linePos = which(abs(sigVals$coeff) == zValue)[1]
	   zVals[linePos] = zValue
	   sigVals$coeff[linePos] = 0
	   remainingLines = which(abs(sigVals$coeff) > eps)
	   selectedLines = c(selectedLines, linePos)
	}

	diffAIC = aic[2:length(aic)] - aic[1:(length(aic)-1)]
	minIndexAIC = which(aic == min(aic))[1]
	nDiffAIC = which(diffAIC > -eps)
	if (length(nDiffAIC) == 0)
	{
	   nDiffAIC = length(aic)
	}else{
	   nDiffAIC = nDiffAIC[1]
	}
	nDiffAIC = nDiffAIC-1
	minIndexAIC = minIndexAIC-1
	nCP = minIndexAIC
	selectedChangePoints = sort(candidateCP2[1:(nCP+2)])
	# for (p in 1:(nCP+1))
	# {
	   # p1 = selectedChangePoints[p]
	   # p2 = selectedChangePoints[p+1]
	# }
	return(selectedChangePoints)
}


####### find change points in y 
findChangePointsForLevel<-function(y, cL)
{
	nL = ceiling(log2(length(y))) # number of levels
	dyadicL = nL -cL # number of levels to use (number of dyadic levels/scales (depth of transform))
	intermidL = 2^((dyadicL-2)+1)-(dyadicL-2)-2 #  number of intermediate scales

	# 1. apply continous haar fitz transform and then shift wavelet coeff. 
	yChft = chft(y,dyadicL,intermidL)
	yChft = shiftWaveletCoeff(yChft,dyadicL,intermidL)

	# 2. find local maximum of wavelet coeff for each scale
	maxima = findWaveletMaximaByScale(yChft)

	# 3. process maxima to get the candidate change points
	threshold = 3
	changePoints = c()
	maximaLines = preprocessMaxima(maxima, yChft, dyadicL ,intermidL, threshold)
	if (length(maximaLines) > 0 & length(ncol(maximaLines)) > 0)
	{
	
		# get most significant scale and position of change 
		sigVals = findSigPosScale(yChft, dyadicL, maximaLines, y)

		# select change points by minimizig AIC
		eps =  2.2204e-16
		changePoints = findOptimizedCP(y, maximaLines, sigVals, eps)
		if (length(changePoints) > 2)
		{
			changePoints = changePoints[2:(length(changePoints)-1)] # remove start and end indices
		}
	}
	return(changePoints)
}

segmentCIM<-function(y, minL = 1, maxL = 1)
{
	changePoints = NULL
	if (minL > maxL | maxL > (0.1*length(y)))
	{
		cat("Check provided minimum and maxiumum levels. The following should hold: minL < maxL < 0.1*length(y)\n")
	}else{
	
		while(minL <= maxL & (length(changePoints) == 0))
		{
			changePoints = findChangePointsForLevel(y, minL)
			minL = minL + 1
		}
	}
	return(list("changePoints" = changePoints, "L" = (minL-1)))
}



