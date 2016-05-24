calcBFMatrix <-function(m1, m2, seg, alpha0 = 1.0)
{
	colnames(seg) = c("chr", "start", "end")
	chrs = as.vector(unique(seg$chr))
	nchrs = length(chrs)
	
	alpha1 = alpha0 + 0.5
	alpha2 = alpha0 + 1
	bfCoeff = (gamma(alpha1)^2)/gamma(alpha2)
	
	#calculate prior mean for cis, trans and self
	cisProp = list("sum" = 0, "n" = 0)
	transProp = list("sum" = 0, "n" =0)
	for (i in 1:nchrs)
	{
		for (j in i:nchrs)
		{
			name1 = chrs[i]
			name2 = chrs[j]
			indicesChr1 = which(seg$chr == name1)
			indicesChr2 = which(seg$chr == name2)
			i1 = indicesChr1[1]
			iN = indicesChr1[length(indicesChr1)]
			j1 = indicesChr2[1]
			jN = indicesChr2[length(indicesChr2)]
			m1Local = m1[i1:iN, j1:jN]
			m2Local = m2[i1:iN, j1:jN]
			
			if (i == j) # cis
			{
				x1 = upperTriangle(m1Local, diag = F)
				x2 = upperTriangle(m2Local, diag = F)
				cisProp$sum = cisProp$sum + sum(x1*2) + sum(x2*2)
				cisProp$n = cisProp$n + length(x1)*4 
			}
			else
			{
				localSum = sum(m1Local + m2Local)
				localN =  (nrow(m1Local)+nrow(m2Local))*(ncol(m1Local)+ncol(m2Local)) 
				transProp$sum = transProp$sum + localSum
				transProp$n = transProp$n + localN
			}
		}
	}  
		
	priorMeanSelf = sum(diag(m1) + diag(m2))/(length(diag(m2))*2)			
	priorMeanCis = cisProp$sum/cisProp$n
	priorMeanTrans = transProp$sum/transProp$n

	#calculate prior variance for cis, trans and self
	sumSquaresCis = 0
	sumSquaresTrans = 0
	for (i in 1:nchrs)
	{
		for (j in i:nchrs)
		{
			name1 = chrs[i]
			name2 = chrs[j]
			indicesChr1 = which(seg$chr == name1)
			indicesChr2 = which(seg$chr == name2)
			i1 = indicesChr1[1]
			iN = indicesChr1[length(indicesChr1)]
			j1 = indicesChr2[1]
			jN = indicesChr2[length(indicesChr2)]
			m1Local = m1[i1:iN, j1:jN]
			m2Local = m2[i1:iN, j1:jN]
			
			if (i == j) # cis
			{
				x1 = upperTriangle(m1Local, diag = F)
				x2 = upperTriangle(m2Local, diag = F)
				localSum = sum((x1 - priorMeanCis)^2)*2 + sum((x2 - priorMeanCis)^2)*2
				sumSquaresCis = sumSquaresCis + localSum
				
			}
			else
			{
				localSum = sum((m1Local - priorMeanTrans)^2) + sum((m2Local - priorMeanTrans)^2)
				sumSquaresTrans = sumSquaresTrans + localSum
			}
		}
	}
	priorVarSelf = (sum((diag(m1) - priorMeanSelf)^2) + sum((diag(m2) - priorMeanSelf)^2))/((length(diag(m2))*2)-1)
	priorVarCis = sumSquaresCis/(cisProp$n-1)
	priorVarTrans =  sumSquaresTrans/(transProp$n-1)

	
	mBfAll = matrix(0, nrow = nrow(m1), ncol = ncol(m1)) # bayes factor matrix
	#compute Bayes Factor for all pairs
	for (i in 1:nchrs)
	{
		for (j in i:nchrs)
		{
			name1 = chrs[i]
			name2 = chrs[j]
			indicesChr1 = which(seg$chr == name1)
			indicesChr2 = which(seg$chr == name2)
			i1 = indicesChr1[1]
			iN = indicesChr1[length(indicesChr1)]
			j1 = indicesChr2[1]
			jN = indicesChr2[length(indicesChr2)]
			m1Local = m1[i1:iN, j1:jN]
			m2Local = m2[i1:iN, j1:jN]
			priorMean = 0
			priorVar = 0
			if (i == j) # cis
			{
				priorMean = priorMeanCis
				priorVar = priorVarCis
			}
			else
			{
				priorMean = priorMeanTrans
				priorVar = priorVarTrans
			}
			obsMean = (m1Local + m2Local)/2
			obsVar = (((m1Local-obsMean)^2) + ((m2Local-obsMean)^2))
			k0 = 2*obsVar/priorVar
			beta0 = 0.5*priorVar
			m1Beta1 = beta0 + k0*((m1Local - priorMean)^2)/(2*k0 + 2)
			m2Beta1 = beta0 + k0*((m2Local - priorMean)^2)/(2*k0 + 2)
			beta2 = beta0 + 0.5*((m1Local - obsMean)^2) + 0.5*((m2Local - obsMean)^2) + (k0*2*((obsMean-priorMean)^2))/(2*k0 + 4)
			mBf = (bfCoeff*(beta0^alpha0)*(beta2^alpha2)*(k0/(k0+1)))/(((k0/(k0+2))^0.5)*(m1Beta1^alpha1)*(m2Beta1^alpha1))
			if (i == j)
			{
				priorMean = priorMeanSelf
				priorVar = priorVarSelf
				obsMean = (m1Local + m2Local)/2
				obsVar = (((m1Local-obsMean)^2) + ((m2Local-obsMean)^2))
				k0 = 2*obsVar/priorVar
				beta0 = 0.5*priorVar
				m1Beta1 = beta0 + k0*((m1Local - priorMean)^2)/(2*k0 + 2)
				m2Beta1 = beta0 + k0*((m2Local - priorMean)^2)/(2*k0 + 2)
				beta2 = beta0 + 0.5*((m1Local - obsMean)^2) + 0.5*((m2Local - obsMean)^2) + (k0*2*((obsMean-priorMean)^2))/(2*k0 + 4)
				m_bf_self = (bfCoeff*(beta0^alpha0)*(beta2^alpha2)*(k0/(k0+1)))/(((k0/(k0+2))^0.5)*(m1Beta1^alpha1)*(m2Beta1^alpha1))
				diag(mBf) = diag(m_bf_self)
			
			}
			mBfAll[i1:iN, j1:jN] = mBf
			
		}
	}
	# note that mBfAll is symmetric so only the upper traingle and the diagonal are filled
	return(mBfAll)

}

compareCIM<-function(m1, m2, seg, bfThreshold = 6.10)
{
	if (nrow(m1) != nrow(seg) | nrow(m1) != nrow(m2) | ncol(m1) != ncol(m2))
	{
		print("provided dimensions of segmentation or matrices do not match")
		return(c())
	}
	
	#standartized version 
	m1S = (m1-mean(m1))/sd(as.numeric(m1))
	m2S = (m2-mean(m2))/sd(as.numeric(m2))
    
	#compute Bayes Factor
	mBf = calcBFMatrix(m1S, m2S, seg[,1:3])
    
	lowerTriangle(mBf) = 0
	chrs = as.vector(unique(seg$chr))
	nchrs = length(chrs)
	
	sigChanges = data.frame()
	sigChangesIndices = which(mBf > bfThreshold, arr.ind = T)
	if (nrow(sigChangesIndices) > 0)
	{ 
		for (k in 1:nrow(sigChangesIndices))
		{
			indexChr1 = as.numeric(sigChangesIndices[k,1])
			indexChr2 = as.numeric(sigChangesIndices[k,2])
			domain1 = seg[indexChr1, ]
			domain2 = seg[indexChr2, ]
			n = length(domain1)
			for (i in 1:n)
			{
				sigChanges[k,i] = domain1[i]	
			}
			for (i in (n+1):(n*2))
			{
				sigChanges[k,i] = domain2[i-n]	
			}
			n = n*2 + 1
			sigChanges[k, n] = mBf[indexChr1, indexChr2]
			sigChanges[k, (n+1)] = m1[indexChr1, indexChr2]
			sigChanges[k, (n+2)] = m2[indexChr1, indexChr2]
			sigChanges[k, (n+3)] = m1S[indexChr1, indexChr2]
			sigChanges[k, (n+4)] = m2S[indexChr1, indexChr2]	
		}			
		# significant changes provides the details (genomic coordinates and any other provided properties) of each pair 
		# along with their Bayes Factor value, and the original and standartized contact frequencies
		colnames(sigChanges) = c(paste(colnames(seg), "1", sep = ""), paste(colnames(seg), "2", sep = ""), "bf", "val1", "val2", "sval1", "sval2")
		sigChanges = sigChanges[order(sigChanges$bf, decreasing = T),]
	}
	# make symmetric
	ind <- lower.tri(mBf)
	mBf[ind] <- t(mBf)[ind]
	return(list("sigChanges" =  sigChanges, "mBf" = mBf))
} 

