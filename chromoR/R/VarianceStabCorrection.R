matrixToVec<-function(m, seg, chrs)
{
	cis = c()
	trans = c()
	iTot = ncol(m)
	nchrs = length(chrs)
	for (i in 1:nchrs)
	{
		indices = which(seg$chr == chrs[i])
		i1 = indices[1]
		iN = indices[length(indices)]
		cis = c(cis, upperTriangle(m[i1:iN, i1:iN]))
		if (iN < iTot)
		{
			trans = c(trans, m[i1:iN, (iN+1):iTot])
		}
	}
	return(list("cis" = cis, "trans" = trans))
}

# padd y with zeros so that it length is a power of 2 
padding <- function(y)
{
	n = length(y)
	if (log2(n) == floor(log2(n)))
	{
		return(y)
	}
	nextPower = floor(log2(n)) + 1
	zeros = rep(0,  ((2^nextPower) - n))
	y = c(y, zeros)
	return(y)
}

# correction procedure - our goal is to have a similar resolution across chromosomes
correctData <- function(y)
{
	# Fisz Haar transform for variance stabilization
	yHft = hft(y)
	# Denoise - shrink coefficients for denoising 
	yHftWd = wd(yHft , filter.number=1, family="DaubExPhase")
	yHftWd = threshold(yHftWd, policy = "cv", dev=madmad)
	yHftWr = wr(yHftWd) # inverse shranked coeff.
	# Apply the inverse Fisz Haar transform to recontrcit the corrected sequence
	return(hft.inv(yHftWr))
}

fillMatrix<-function(m, cis, trans, seg, chrs)
{
	nchrs = length(chrs)
	iTot = ncol(m)
	lastCisIndex = 1
	lastTransIndex = 1
	for (i in 1:nchrs)
	{
		indices = which(seg$chr == chrs[i])
		i1 = indices[1]
		iN = indices[length(indices)]
		nDomains = iN-i1+1
		size = nDomains^2
		
		# fill cis
		currCisSize = ((size) - sqrt(size))/2
		s1 = lastCisIndex
		sN = lastCisIndex + currCisSize -1 
		vCis = cis[s1:sN]
		lastCisIndex = sN + 1
		mCis = matrix(0, nrow = nDomains, ncol = nDomains)
		upperTriangle(mCis) = vCis
		m[i1:iN, i1:iN] = mCis
		
		# fill trans
		if (iN < iTot)
		{
			currTransSize = (iTot-iN)*nDomains
			s1 = lastTransIndex
			sN = lastTransIndex + currTransSize -1 
			vTrans = trans[s1:sN]
			lastTransIndex = sN + 1
			mTrans = matrix(vTrans)
			m[i1:iN, (iN+1):iTot] = mTrans
		}
	}
	return(m)
}

correctCIM<-function(m, seg, removeUncovered = FALSE)
{
	if (nrow(m) != nrow(seg))
	{
		print("provided segmentation does not match contact matrix")
		return(c())
	}
	colnames(seg) = c("chr", "start", "end")
	chrs = as.vector(unique(seg$chr))
	nchrs = length(chrs)
	
	if (removeUncovered)
	{
		x = rowSums(m) 
		ind = which(x == 0)
		seg = cbind(seg, x)
		seg = subset(seg, seg$x > 0)
		m = m[!(apply(m, 1, function(y) all(y == 0))),]
		m = m[,!(apply(m, 2, function(y) all(y == 0)))]
		seg = seg[,1:3]
	}
	
	
	# convert to trans, cis and self vectors
	pairwise = matrixToVec(m, seg, chrs)
	cis = pairwise$cis
	trans = pairwise$trans
	self = diag(m)

	# correction procedure 
	# self
	self = padding(self)
	correctedSelf = correctData(self)

	# cis: apply correction for each pair separatedly 
	lastCisIndex = 1
	correctedCis = c()
	for (i in 1:nchrs)
	{
		indices = which(seg$chr == chrs[i])
		i1 = indices[1]
		iN = indices[length(indices)] 
		mSize = (iN-i1+1)^2
		# fill cis
		vSize = ((mSize) - sqrt(mSize))/2
		s1 = lastCisIndex
		sN = lastCisIndex + vSize -1 
		currCis = cis[s1:sN]
		lastCisIndex = sN + 1
		origSize = length(currCis)
		currCis = padding(currCis)
		currCorrectedCis = correctData(currCis)
		correctedCis = c(correctedCis, currCorrectedCis[1:origSize]) 
	}

	# trans: apply correction for each pair separatedly
	correctedTrans = c()
	lastTransIndex1 = 1
	iTot = nrow(m)
	for (i in 1:nchrs)
	{
		indices = which(seg$chr == chrs[i])
		i1 = indices[1]
		iN = indices[length(indices)]
		nRegions = iN-i1+1
		if (iN < iTot)
		{
			vSize = (iTot-iN)*nRegions
			s1 = lastTransIndex1
			sN = lastTransIndex1 + vSize -1 
			currTrans = trans[s1:sN]
			currCorrectedTrans = c()
			k = i+1
			lastTransIndex2 = 1
			for (j in k:nchrs)
			{
				indices = which(seg$chr == chrs[j])
				j1 = lastTransIndex2
				jN = j1 + length(indices)*nRegions
				if (jN > length(currTrans))
				{
					jN = length(currTrans)
				}
				
				currPairwiseTrans = currTrans[j1:jN]
				origSize = length(currPairwiseTrans)
				currPairwiseTrans = padding(currPairwiseTrans)	
				correctedPairwiseTrans =  correctData(currPairwiseTrans)	
				currCorrectedTrans =  c(currCorrectedTrans, correctedPairwiseTrans[1:origSize])	
				lastTransIndex2 = jN + 1	
			}
			lastTransIndex1 = sN + 1
			correctedTrans = c(correctedTrans, currCorrectedTrans)
		}
	}

	# reconstruct the matrix from the sequences
	mCorrected = matrix(0, ncol = ncol(m), nrow = nrow(m))
	mCorrected = fillMatrix(mCorrected, correctedCis, correctedTrans, seg, chrs)

	# make symmetric with zeros in diagonal
	ind <- lower.tri(mCorrected)
	mCorrected[ind] <- t(mCorrected)[ind]
	diag(mCorrected) = correctedSelf[1:length(diag(m))]
	
	return(list("mCorrected" = mCorrected, "seg" = seg))
		
}

correctPairCIM<-function(m, isCis)
{
	mCorrected = matrix(0, nrow = nrow(m), ncol = ncol(m))
	if (isCis)
	{
		self = diag(m)
		self = padding(self)
		correctedSelf = correctData(self)
		
		cis = upperTriangle(m, diag = F)
		origSize = length(cis)
		cis = padding(cis)
		correctedCis = correctData(cis)
		
		
		upperTriangle(mCorrected) = correctedCis[1:origSize]
		ind <- lower.tri(mCorrected)
		mCorrected[ind] <- t(mCorrected)[ind]
		diag(mCorrected) = correctedSelf[1:length(diag(m))]
	
	}else{
		trans = as.vector(m, mode = "numeric")
		origSize = length(trans)
		trans = padding(trans)
		correctedTrans = correctData(trans)
		correctedTrans = correctedTrans[1:origSize]	
		mCorrected = matrix(correctedTrans, nrow = nrow(m), ncol = ncol(m))
	}
	return(mCorrected)
}














