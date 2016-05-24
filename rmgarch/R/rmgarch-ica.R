#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2013.
##   This file is part of the R package rmgarch.
##
##   The R package rmgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rmgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

#------------------------------------------------------------------------------
# ICA Algorithms/GO-GARCH
#------------------------------------------------------------------------------
.makeiid = function(X, method = "fastica", ica.fix = list(), demean = FALSE, ...)
{
	X 	= as.matrix(X)
	n 	= NROW(X)
	m 	= NCOL(X)
	# should modify to deal with non square matrix where they are all passed in list
	if(!is.null(ica.fix$A)){
		A = as.matrix(ica.fix$A)
		# reduced dimension model
		if(NROW(A) != NCOL(A)){
			fc = NCOL(A)
			if(fc > NROW(A)) stop("\nA matrix provided in ica.fix must be no.signals(m) by n.factors(f) where f<=m.")
			if(!is.null(ica.fix$W)){
				W = as.matrix(ica.fix$W)
				if(NROW(W) != fc | NCOL(W) != m) stop("\nwrong dimension for W matrix provided in ica.fix\n", call.  = FALSE)
			} else{
				stop("\nW matrix must be provided in reduced dimension system.")
			}
			if(!is.null(ica.fix$K)){
				K = as.matrix(ica.fix$K)
				if(NROW(K) != fc | NCOL(K) != m) stop("\nwrong dimension for K matrix provided in ica.fix\n", call.  = FALSE)
			} else{
				stop("\nK (whiteningMatrix) matrix must be provided in reduced dimension system.")
			}
			if(!is.null(ica.fix$Kinv)){
				Kinv = as.matrix(ica.fix$Kinv)
				if(NROW(Kinv) != m | NCOL(Kinv) != fc) stop("\nwrong dimension for Kinv matrix provided in ica.fix\n", call.  = FALSE)
			} else{
				stop("\nKinv (dewhiteningMatrix) matrix must be provided in reduced dimension system.")
			}
			U =  t( W %*% Kinv )
			# equivalent: t( t(A) %*% t(k) )
			if(demean){
				X = scale(X, scale = FALSE)
				mu = attr(X, "scaled:center")
			} else{
				mu = rep(0, m)
			}
			Y = X %*% t(W) + (matrix(1, nrow = NROW(X), ncol = 1)%*% (mu %*% t(W)))
		} else{
			if(NROW(A) != m | NCOL(A) != m) stop("\nwrong dimension for A matrix provided in ica.fix\n", call.  = FALSE)
			if(!is.null(ica.fix$K)){
				K = as.matrix(ica.fix$K)
				if(NROW(K) != m | NCOL(K) != m) stop("\nwrong dimension for K matrix provided in ica.fix\n", call.  = FALSE)
			} else{
				stop("\nK (whiteningMatrix) matrix must be provided")
			}
			# ToDo: checks for W
			if(is.null(ica.fix$W)) W = solve(A) else W = ica.fix$W
			K = K
			Kinv = solve(K)
			if(demean){
				X = scale(X, scale = FALSE)
				mu = attr(X, "scaled:center")
			} else{
				mu = rep(0, m)
			}
			Y = X %*% t(W) + (matrix(1, nrow = NROW(X), ncol = 1)%*% (mu %*% t(W)))
			U =  t( W %*% Kinv )
		}
		return(list(Y = Y, W = W, A = A, U = U, K = K, Kinv = Kinv, E = NULL, D = NULL))
	} else{
		ica 	= switch(method,
				fastica = .fastica(X, demean = demean, ...),
				radical = .radical(X, demean = demean, ...))
		W = ica$W
		K = ica$K
		A = ica$A
		Y = ica$Y 
		Kinv = ica$Kinv
		U = ica$U
		D = ica$D
		E = ica$E
		return(list(Y = Y, W = W, A = A, U = U, K = K, Kinv = Kinv, E = E, D = D))
	}
}


.fastica = function(x, ...)
{	
	ans = fastica(X = x, ...)
	W = ans$W
	A = ans$A
	Y = ans$S
	K = ans$whiteningMatrix
	Kinv = ans$dewhiteningMatrix
	E = ans$E
	D = ans$D
	U = ans$U
	return(list(W = W, A = A, U = U, Y = Y, K = K, Kinv = Kinv, E = E, D = D))
}

#.fastICA = function(x, ...)
#{	
#	ans = fastICA(X = x, n.comp = ncol(x), method = "C", ...)	
#	W = ans$K %*% ans$W
#	A = ans$A
#	Y = ans$S
#	K = ans$whiteningMatrix
#	Kinv = ans$dewhiteningMatrix
#	E = ans$E
#	D = ans$D
#	U = ans$U
#	return(list(W = W, A = A, U = U, Y = Y, K = K, Kinv = Kinv, E = E, D = D))
#}
.radical = function(x, ...)
{
	ans = radical(X = x, ...)
	W = ans$W
	A = ans$A
	Y = ans$S
	K = ans$whiteningMatrix
	Kinv = ans$dewhiteningMatrix
	E = ans$E
	D = ans$D
	U = ans$U
	return(list(W = W, A = A, U = U, Y = Y, K = K, Kinv = Kinv, E = E, D = D))
}
#############################################################################
# FASTICA - Fast Independent Component Analysis
# FastICA for Matlab 7.x and 6.x
# Version 2.5, October 19 2005
# Copyright (c) Hugo Gavert, Jarmo Hurri, Jaakko Sarela, and Aapo Hyvarinen.
# A translation of the matlab function from the FastICA toolbox
# http://www.cis.hut.fi/projects/ica/fastica/
# R translation: 2009 Alexios Ghalanos
# 			   : 2012 Extra Methods for pca.cov
#############################################################################

fastica = function(X, approach = c("symmetric", "deflation"), n.comp = dim(X)[2], demean = TRUE,
		pca.cov = c("ML", "LW", "ROB", "EWMA"), gfun = c("pow3", "tanh", "gauss", "skew"), 
		finetune = c("none", "pow3", "tanh", "gauss", "skew"),
		tanh.par = 1, gauss.par = 1, step.size = 1, stabilization = FALSE, 
		epsilon = 1e-4, maxiter1 = 1000, maxiter2 = 5, A.init = NULL, pct.sample = 1, 
		firstEig = NULL, lastEig = NULL, pcaE = NULL, pcaD = NULL, 
		whiteSig = NULL, whiteMat = NULL, dewhiteMat = NULL, rseed = NULL, trace = FALSE, ...)
{
	########################################
	tic = Sys.time()
	# force evaluation of n.comp before transposing else will lazily evaluate
	# afterwards and we'll get the wrong result
	n.comp = n.comp
	# transpose dataset (as expected by most ICA algorithms)
	X = t(X)
	idx = match(tolower(pca.cov[1]), tolower(c("ML", "LW", "ROB", "EWMA")))
	if(is.na(idx)) stop("\nfastica-->error: pca.cov choice not recognized.")
	pca.cov = c("ML", "LW", "ROB", "EWMA")[idx]
	
	########################################
	# STAGE 1: Demean Data
	
	# Remove the mean and check the data
	if(demean){ 
		mixedmean = apply(X, 1, "mean")
		mixedsig  = t(apply(X, 1, FUN = function(x) x - mean(x) ))
	} else{
		mixedmean = rep(0, dim(X)[1])
		mixedsig  = X
	}
	Dim = dim(mixedsig)[1]
	nsamples = dim(mixedsig)[2]
	myy = step.size
	
	
	# n.comps + firstEig/lastEig choice
	if(n.comp>nrow(X)) stop("\nCannot Choose more components than signals!")
	if(is.null(firstEig)) firstEig = 1
	if(is.null(lastEig)) lastEig = n.comp
	if(!is.null(firstEig)){
		if(is.null(lastEig)){
			lastEig = firstEig+n.comp
			if(lastEig>nrow(X)) stop("\nfirstEig and n.comp combination greater than no. of signals.")
		} else{
			if(lastEig>nrow(X)) stop("\nlastEig greater than no. of signals.")
			n.comp = lastEig - firstEig + 1
			if(n.comp<1) stop("\nCannot have zero factors in PCA!")
		}
	}
	if(!is.null(lastEig)){
		n.comp = lastEig - firstEig + 1
	}
	# PCA E (m x f)
	if(!is.null(pcaE)){
		if(is.null(pcaD)) stop("\npcaD required when passing pcaE.")
		if(!is.matrix(pcaE)) stop("\npcaE must be an m by f matrix (see documentation)")
		if(nrow(pcaE)!=nrow(X)) stop("\npcaE no.rows must equal no.signals.")
		if(ncol(pcaE)!=n.comp) stop("\npcaE no.cols must equal no.comps.")
	}
	# PCA D (f by f)
	if(!is.null(pcaD)){
		if(is.null(pcaE)) stop("\npcaE required when passing pcaD.")
		if(!is.matrix(pcaD)) stop("\npcaD must be an f by f diagonal matrix (see documentation)")
		if(nrow(pcaD)!=n.comp) stop("\npcaD no.rows must equal no.comps.")
		if(ncol(pcaD)!=n.comp) stop("\npcaD no.cols must equal no.comps.")
		# recover the covariance matrix (for non-dimensionality reduced systems)
		if(NROW(pcaE) == NCOL(pcaE)) C = pcaE %*% pcaD %*% solve(pcaE) else C = NULL
	}
	
	# whiteSig (n by f)
	if(!is.null(whiteSig)){
		if(!is.matrix(whiteSig)) stop("\nwhiteSig must be an n by f matrix (see documentation)")
		if(nrow(whiteSig)!=ncol(X)) stop("\nwhiteSig no.rows must equal no.rows Data.")
		if(ncol(whiteSig)!=n.comp) stop("\nwhiteSig no.cols must equal no.comps.")
		# now transpose
	} else{
		if(!is.null(pcaE)){
			whiteMat = solve(sqrt(pcaD))%*%t(pcaE)
			whiteSig = t(mixedsig) %*% t(whiteMat)
			whiteSig = t(whiteSig)
		}
	}
	
	if(!is.null(dewhiteMat)){
		if(is.null(whiteMat)) stop("\ndewhiteMat cannot be provided without whiteMat.")
		if(nrow(dewhiteMat)!=nrow(X)) stop("\ndewhiteMat must be an m by f matrix (see documentation)")
		if(ncol(dewhiteMat)!=n.comp) stop("\ndewhiteMat must be an m by f matrix (see documentation)")
	}
	
	# whiteMat (m by f)
	if(!is.null(whiteMat)){
		if(!is.matrix(whiteMat)) stop("\nwhiteMat must be an f by m matrix (see documentation)")
		if(nrow(whiteMat)!=n.comp) stop("\nwhiteMat no.rows must equal no.comps.")
		if(ncol(whiteMat)!=nrow(X)) stop("\nwhiteMat no.cols must equal no.signals.")
		if(is.null(dewhiteMat)){
			if(ncol(whiteMat)==nrow(whiteMat)){
				dewhiteMat = t(solve(whiteMat))
			} else{
				stop("\ndewhiteMat needs to be supplied for user supplied non-symmetric whiteMat")
			}
		}
		# now transpose
		if(is.null(whiteSig)){
			whiteSig = t(mixedsig) %*% t(whiteMat)
			whiteSig = t(whiteSig)
		}
	}
	
	
	# some basic checks
	if(length(which(is.na(X)))>0) stop("\nfastica-->error: NA's in data. Remove and retry.\n", call. = FALSE)
	
	approach = approach[1]
	if(!any(approach == c("symmetric", "deflation"))) stop("\nfastica-->error: Invalid input for approach.\n", call. = FALSE)
	
	gfun = tolower(gfun)[1]
	if(!any(gfun == c("pow3", "tanh", "gauss", "skew"))) stop("\nfastica-->error: Invalid input for gfun.\n", call. = FALSE)
	
	finetune = tolower(finetune)[1]
	if(!any(finetune == c("pow3", "tanh", "gauss", "skew", "none"))) stop("\nfastica-->error: Invalid input for finetune.\n", call. = FALSE)
	
	########################################
	# STAGE 2: PCA decomposition to obtain eigenvalues and eigenvectors, 
	# and optionally reduce dimension
	
	# PCA
	jumpPCA = 0
	E = NULL
	D = NULL
	whiteningsig = NULL
	whiteningMatrix = NULL
	dewhiteningMatrix = NULL
	if(!is.null(pcaE))
	{
		jumpPCA = jumpPCA + 1;
		E = pcaE
	}
	if(!is.null(pcaD))
	{
		jumpPCA = jumpPCA + 1;
		D = pcaD
	}
	# calculate if there are enought parameters to skip PCA and whitening
	jumpWhitening = 0
	if(!is.null(whiteSig))
	{
		jumpWhitening = jumpWhitening + 1
		whiteningsig = whiteSig
	}
	if(!is.null(whiteMat))
	{
		jumpWhitening = jumpWhitening + 1
		whiteningMatrix = whiteMat
	}
	if(!is.null(dewhiteMat))
	{
		jumpWhitening = jumpWhitening + 1
		dewhiteningMatrix = dewhiteMat
	}
	if(trace){
		cat(paste("\nNo. of signals: ", Dim, sep = ""))
		cat(paste("\nNo. of samples: ", nsamples, "\n", sep = ""))
	}
	if(Dim > nsamples && trace) warning("\nfastica-->warning: The signal matrix(X) may be orientated in the wrong way...\n")
	if(jumpWhitening == 3){
		if(trace) cat("\nWhitened signal and corresponding matrices supplied. PCA calculations not needed.\n")
	} else{
		if(jumpPCA == 2) {
			if(trace) cat("\nValues for PCA calculations supplied. PCA calculations not needed.\n")
		} else{
			if(jumpPCA == 1 && trace) cat("\nYou must supply both pcaE and pcaD in order to jump PCA.\n")
			tmp = .pca(mixedsig, firstEig, lastEig, pca.cov, trace, ...)
			# D = diag(eigenvalues)
			# E = tmp$vectors
			# C = covariance matrix
			E = tmp$E
			D = tmp$D
			C = tmp$C
		}
	}
	########################################
	# STAGE 3: Use PCA values to whiten/decorrelate the Data and return the 
	# whitening and de-whitening matrix
	# 	whiteningMatrix = solve(sqrt(D))%*%t(E)
	#   dewhiteningMatrix = E %*% sqrt(D)
	
	if(jumpWhitening == 3){
		if(trace) cat("\nWhitening not needed.\n")
	} else{
		tmp = whitenv(mixedsig, E, D, trace)
		whiteningsig = tmp$newVectors
		whiteningMatrix = tmp$whiteningMatrix
		dewhiteningMatrix = tmp$dewhiteningMatrix
		#print(whiteningMatrix, digits = 5)
	}
	
	Dim = dim(whiteningsig)[1]
	
	if(n.comp > Dim) 
	{
		n.comp = Dim
		# Show warning only if verbose = 'on' and user supplied a value for 'n.comp'
		if (trace){
			warning(paste("\nfastica-->warning: estimating only ", n.comp, " independent components\n", sep = "" ))
			cat("(Can't estimate more independent components than dimension of data)\n")
		}
	}
	########################################
	# STAGE 4: ICA rotation to seperate out independent signals. Return mixing (A) and unmixing matrices (W)
	# Calculate the ICA with fixed point algorithm.
	tmp = .fpica(whiteningsig, whiteningMatrix, dewhiteningMatrix, approach, n.comp, gfun, finetune, tanh.par, gauss.par, myy, 
			stabilization, epsilon, maxiter1, maxiter2, A.init, pct.sample, rseed, trace)
	A = tmp$A
	W = tmp$W
	# Generate Signals (generative model)
	if(!is.null(A)){
		S = W %*% mixedsig + (W %*% mixedmean) %*% .ones(1, nsamples)
	} else{
		S = NULL
	}
	U = t(A) %*% t(tmp$whiteningMatrix)
	
	toc = Sys.time() - tic
	
	return(list(A = A, W = W, U = t(U), S = t(S), E = E, D = D, C = C, mu = mixedmean, 
					whiteningMatrix = (tmp$whiteningMatrix), 
					dewhiteningMatrix = (tmp$dewhiteningMatrix), 
					rseed = tmp$rseed, elapsed = toc))
}
# A translation of the matlab function from the FastICA toolbox
# http://www.cis.hut.fi/projects/ica/fastica/

.fpica = function(X, whiteningMatrix, dewhiteningMatrix, approach = c("symmetric","deflation"), 
		n.comp = dim(X)[2], gfun = c("pow3", "tanh", "gauss", "skew"), finetune = c("none", "pow3", "tanh", "gauss", "skew"), 
		tanh.par = 1, gauss.par = 1, myy = 1, stabilization = TRUE, epsilon = 1e-4, maxiter1 = 1000, 
		maxiter2 = 100, A.init = NULL, pct.sample = 1, rseed = NULL, trace = FALSE)
{
	vectorSize = dim(X)[1]
	nsamples =  dim(X)[2]
	# check inputs
	if(is.null(rseed)) rseed = runif(1, 1, 1000000)
	# n.comp
	if(vectorSize < n.comp)
		stop("\nfastica-->error: n.comp is greater than dimension of X!\n", call. = FALSE)
	# approach
	approach = approach[1]
	if(!any(tolower(approach) == c("symmetric","deflation")))
		stop("\nfastica-->error: invalid input value for approach.\n", call. = FALSE)
	
	# gfun
	if(!any(tolower(gfun) == c("pow3", "tanh", "gauss", "skew")))
		stop("\nfastica-->error: invalid input value for gfun.\n", call. = FALSE)
	# sampleSize
	if(pct.sample > 1) {
		pct.sample = 1
		if(trace) warning("\nfpica-->warning: Setting pct.sample to 1.\n")
	} else{
		if ( (pct.sample * nsamples) < 1000) {
			pct.sample = min(1000/nsamples, 1)
			if(trace) warning(paste("\nfastica-->warning: Setting pct.sample to", round(pct.sample,2), "(",
								floor(pct.sample * nsamples)," % of sample).\n"), sep = "")
		}
	}
	# Checking the value for nonlinearity.
	gOrig = switch(tolower(gfun),
			pow3  = 10,
			tanh  = 20,
			gauss = 30,
			skew  = 40)
	if(pct.sample != 1) gOrig = gOrig + 2
	if(myy != 1) gOrig = gOrig + 1
	
	finetune = finetune[1]
	finetuningEnabled = 1
	gFine = switch(tolower(finetune),
			pow3  = 10 + 1,
			tanh  = 20 + 1,
			gauss = 30 + 1,
			skew  = 40 + 1,
			none  = if(myy != 1) gOrig else gOrig + 1)
	
	if(tolower(finetune) == "none") finetuningEnabled = 0
	
	if(stabilization){
		stabilizationEnabled = 1
	} else{
		if(myy != 1) stabilizationEnabled = 1 else stabilizationEnabled = 0
	}
	
	myyOrig = myy
	#% When we start fine-tuning we'll set myy = myyK * myy
	myyK = 0.01
	# How many times do we try for convergence until we give up.
	failureLimit = 5
	usedNlinearity = gOrig
	stroke = 0
	notFine = 1
	long = 0
	
	# Checking the value for initial state.
	if(is.null(A.init)){
		initialStateMode = 0
		if(trace) cat("\nfastica-->Using random initialization.\n")
	} else{
		if(dim(A.init)[1] != dim(whiteningMatrix)[2]){
			initialStateMode = 0
			A.init = NULL
			if(trace) warning("\nfastica-->warning: size of initial guess is incorrect. Using random initial guess.\n")
		} else{
			initialStateMode = 1
			if(dim(A.init)[2] < n.comp){
				if(trace) warning(paste("\nfastica-->warning: initial guess only for first ", dim(A.init)[2]," components. Using random initial A.init for others.\n", sep = ""))
				set.seed(rseed)
				A.init[,dim(A.init)[2] + (1:n.comp)]  = matrix(runif(vectorSize * (n.comp-dim(A.init)[2]))-.5, ncol = n.comp-dim(A.init)[2])
			} else{
				if(dim(A.init)[2] > n.comp && trace) warning(paste("\nfastica-->warning: Initial guess too large. The excess column are dropped.\n", sep = ""))
				A.init = A.init[, 1:n.comp]				
			}
			if(trace) cat("\nfastica-->Using initial A.init.\n")
		}
	}
	
	if(trace) cat("\nStarting ICA calculation...\n")
	
	# symmetric case
	ans = switch(approach,
			symmetric = .fpsymm(X, gOrig, n.comp, vectorSize, initialStateMode, whiteningMatrix, dewhiteningMatrix, A.init, 
					maxiter1, epsilon, stabilizationEnabled, finetuningEnabled, failureLimit, myy, myyK, myyOrig, nsamples, 
					pct.sample, tanh.par, gauss.par, gFine, trace, rseed),
			deflation = .fpdefl(X, gOrig, n.comp, vectorSize, initialStateMode, whiteningMatrix, dewhiteningMatrix, A.init, 
					maxiter1, maxiter2, epsilon, stabilizationEnabled, finetuningEnabled, failureLimit, myy, myyK, myyOrig, 
					nsamples, pct.sample, tanh.par, gauss.par, gFine, trace, rseed))
	
	return(list(A = ans$A, W = ans$W, whiteningMatrix = whiteningMatrix, dewhiteningMatrix = dewhiteningMatrix,
					rseed = rseed))
}

.fpsymm = function(X, gOrig, n.comp, vectorSize, initialStateMode, whiteningMatrix, dewhiteningMatrix, A.init, 
		maxiter1, epsilon, stabilizationEnabled, finetuningEnabled, failureLimit, myy, myyK, myyOrig, nsamples, pct.sample, 
		tanh.par, gauss.par, gFine, trace, rseed)
{
	usedNlinearity = gOrig
	stroke = 0
	notFine = 1
	long = 0
	A = .zeros(vectorSize, n.comp)  # Dewhitened basis vectors.
	if(initialStateMode == 0){
		# Take random orthonormal initial vectors.
		set.seed(rseed)
		B = .orthogonal(matrix(rnorm(vectorSize*n.comp), ncol = n.comp, byrow = TRUE))
	} else{
		# Use the given initial vector as the initial state
		B = whiteningMatrix %*% A.init
	}
	BOld2 = BOld  = .zeros(dim(B)[1], dim(B)[2])
	# This is the actual fixed-point iteration loop.
	for(i in 1:(maxiter1 + 1))
	{
		if(i == maxiter1 + 1){
			cat(paste("No convergence after ", maxiter1," steps\n", sep = ""))
			if(!is.null(B)){
				B = B * .sqrtm(Re(solve(t(B) %*% B)))
				W = t(B) %*% whiteningMatrix
				A = dewhiteningMatrix %*% B
			} else{
				W = NULL
				A = NULL
			}
			# exit condition
			return(list(A = A, W = W))
		}
		# Symmetric orthogonalization.
		B = B %*% Re(.sqrtm(solve(t(B) %*% B)))
		# Test for termination condition. Note that we consider opposite
		# directions here as well.
		minAbsCos = min(abs(diag(t(B) %*% BOld)))
		minAbsCos2 = min(abs(diag(t(B) %*% BOld2)))
		if( (1 - minAbsCos) < epsilon ){
			if(finetuningEnabled && notFine){
				if(trace) cat("\nInitial convergence, fine-tuning: \n")
				notFine = 0
				usedNlinearity = gFine
				myy = myyK * myyOrig
				BOld  = .zeros(dim(B)[1], dim(B)[2])
				BOld2 = .zeros(dim(B)[1], dim(B)[2])
			} else{
				if(trace) cat(paste("\nConvergence after ", i, "steps\n", sep = ""))
				# Calculate the de-whitened vectors.
				A = dewhiteningMatrix %*% B
				break()
			}
			
		}
		if(stabilizationEnabled)
		{
			if(!stroke)
			{
				if( (1 - minAbsCos2) < epsilon )
				{
					if(trace) cat("\nStroke!\n")
					stroke = myy
					myy = .5*myy;
					if( usedNlinearity%%2 == 0) usedNlinearity = usedNlinearity + 1
				}
			} else{
				myy = stroke
				stroke = 0
				if( (myy == 1) && (usedNlinearity%%2 != 0) ) usedNlinearity = usedNlinearity - 1
			}
			if(!long && (i>maxiter1/2) )
			{
				if(trace) cat("\nTaking long (reducing step size)\n")
				long = 1
				myy = .5 * myy
				if(usedNlinearity%%2 == 0) usedNlinearity = usedNlinearity + 1
			}
		}
		BOld2 = BOld
		BOld = B
		# Show the progress...
		if(trace)
		{
			if(i == 1)
			{
				cat(paste("Step no. ", i,"\n", sep = ""))
			} else{
				cat(paste("Step no. ", i,", change in value of estimate ", round(1 - minAbsCos,3),"\n", sep = ""))				
			}
		}
		#A = dewhiteningMatrix %*% B
		#W = t(B) %*% whiteningMatrix
		# Nlinearity cases
		fmethod = paste("f", usedNlinearity, sep = "")
		B = .fsmethod(fmethod, X, B, myy, nsamples, pct.sample, tanh.par, gauss.par)
	}
	# Calculate ICA filters.
	W = t(B) %*% whiteningMatrix
	return(list(W = W, A = A))
}

.fpdefl = function(X, gOrig, n.comp, vectorSize, initialStateMode, whiteningMatrix, dewhiteningMatrix, A.init, 
		maxiter1, maxiter2, epsilon, stabilizationEnabled, finetuningEnabled, failureLimit, myy, myyK, myyOrig, nsamples, pct.sample, 
		tanh.par, gauss.par, gFine, trace, rseed)
{
	B = .zeros(vectorSize, vectorSize)
	A = .zeros(dim(dewhiteningMatrix)[1], dim(dewhiteningMatrix)[2])
	W = .zeros(dim(dewhiteningMatrix)[2], dim(dewhiteningMatrix)[1])
	# The search for a basis vector is repeated n.comp times.
	j = 1
	numFailures = 0
	while(j <= n.comp){
		myy = myyOrig
		usedNlinearity = gOrig
		stroke = 0
		notFine = 1
		long = 0
		endFinetuning = 0
		# Show the progress...
		if(trace) cat(paste("IC", j, sep =""))
		# Take a random initial vector of lenght 1 and orthogonalize it
		# with respect to the other vectors.
		if(is.null(A.init)){
			set.seed(rseed)
			w = matrix(rnorm(vectorSize), ncol = 1)
		} else{
			w = whiteningMatrix %*% A.init[,j]
		}
		w = w - B %*% t(B) %*% w
		w = w / .fnorm(w)
		wOld2 = wOld = .zeros(dim(w)[1], dim(w)[2])
		# This is the actual fixed-point iteration loop.
		#    for i = 1 : maxNumIterations + 1
		i = 1
		gabba = 1
		while(i <= (maxiter1 + gabba)){
			# Project the vector into the space orthogonal to the space
			# spanned by the earlier found basis vectors. Note that we can do
			# the projection with matrix B, since the zero entries do not
			# contribute to the projection.
			w = w - B %*% t(B) %*% w
			w = w / .fnorm(w)
			if(notFine){
				if(i == (maxiter1 + 1)){
					if(trace) cat(paste("\nComponent number ", j, " did not converge in ", maxiter1, " iterations.\n", sep = ""))
					j = j - 1
					numFailures = numFailures + 1
					if(numFailures > failureLimit){
						if(trace)  cat(paste("\nfastica-->: error: Too many failures to converge (",  numFailures, "). Giving up.\n", sep = ""))
						if(j == 0){
							A = NULL
							W = NULL
						}
						return(list(W = W, A = A))
					}
				} 
			} else{
				if(i >= endFinetuning) wOld = w
			}
			if(trace) cat(".")
			# Test for termination condition. Note that the algorithm has
			# converged if the direction of w and wOld is the same, this
			# is why we test the two cases.
			if( (.fnorm(w - wOld) < epsilon) | (.fnorm(w + wOld) < epsilon) ){
				if(finetuningEnabled && notFine){
					if(trace) cat("\nInitial convergence, fine-tuning\n")
					notFine = 0
					gabba = maxiter2
					wOld2 = wOld = .zeros(dim(w)[1], dim(w)[2])
					usedNlinearity = gFine
					myy = myyK * myyOrig
					endFinetuning = maxiter2 + i
				} else{
					numFailures = 0
					# Save the vector
					B[ , j] = w
					# Calculate the de-whitened vector.
					A[ , j] = dewhiteningMatrix %*% w
					# Calculate ICA filter.
					W[j,] = t(w) %*% whiteningMatrix
					if(trace) cat(paste("\ncomputed ( ", i, " steps ) \n", sep = ""))
					break()
				}
			}
			if(stabilizationEnabled){
				if(!stroke && ( (.fnorm(w - wOld2) < epsilon) | (.fnorm(w + wOld2) < epsilon) )){
					stroke = myy
					if(trace) cat("\nStroke!\n")
					myy = .5 * myy
					if( (usedNlinearity%%2) == 0) usedNlinearity = usedNlinearity + 1
				}
				if(stroke){
					myy = stroke
					stroke = 0
					if( (myy == 1) && ( (usedNlinearity%%2) != 0)) usedNlinearity = usedNlinearity - 1
				}
				if(notFine && !long && (i > maxiter1/2)){
					if(trace) cat("\nTaking too long (reducing step size\n")
					long = 1
					myy = .5 * myy
					if( (usedNlinearity%%2) == 0) usedNlinearity = usedNlinearity + 1
				}	
			}
			wOld2 = wOld
			wOld = w
			fmethod = paste("f", usedNlinearity, sep = "")
			#print(fmethod)
			w = .fdmethod(fmethod, X, w, myy, nsamples, pct.sample, tanh.par, gauss.par)
			w = w / .fnorm(w)
			i = i + 1
		}
		j = j + 1
	}
	if(trace) cat("\nDone.\n")
	return(list(W = W, A = A))
}

#--------------------------------------------------------------------------------------------------------
# symmetric methods
#--------------------------------------------------------------------------------------------------------
.fsmethod = function(fmethod, X, B, myy, nsamples, pct.sample, tanh.par, gauss.par)
{
	ans = switch(fmethod,
			f10 = .fs10(X, B,  myy, nsamples,  pct.sample),
			f11 = .fs11(X, B,  myy, nsamples,  pct.sample),
			f12 = .fs12(X, B,  myy, nsamples,  pct.sample),
			f13 = .fs13(X, B,  myy, nsamples,  pct.sample),
			f20 = .fs20(X, B,  myy, tanh.par,  nsamples, pct.sample),
			f21 = .fs21(X, B,  myy, tanh.par,  nsamples, pct.sample),
			f22 = .fs22(X, B,  myy, tanh.par,  nsamples, pct.sample),
			f23 = .fs23(X, B,  myy, tanh.par,  nsamples, pct.sample),
			f30 = .fs30(X, B,  myy, gauss.par, nsamples, pct.sample),
			f31 = .fs31(X, B,  myy, gauss.par, nsamples, pct.sample),
			f32 = .fs32(X, B,  myy, gauss.par, nsamples, pct.sample),
			f33 = .fs33(X, B,  myy, gauss.par, nsamples, pct.sample),
			f40 = .fs40(X, B,  myy, nsamples,  pct.sample),
			f41 = .fs41(X, B,  myy, nsamples,  pct.sample),
			f42 = .fs42(X, B,  myy, nsamples,  pct.sample),
			f43 = .fs43(X, B,  myy, nsamples,  pct.sample))
	return(ans)
}

# pow3
.fs10 = function(X, B, myy, nsamples, pct.sample)
{
	ans = (X %*% (( t(X) %*% B)^ 3)) / nsamples - 3 * B
	ans
}

.fs11 = function(X, B, myy, nsamples, pct.sample)
{
	Y = t(X) %*% B
	Gpow3 = Y^3
	Beta = as.numeric(apply(Y * Gpow3, 2, "sum"))
	D = diag(1 / (Beta - 3 * nsamples))
	ans = B + myy * B %*% (t(Y) %*% Gpow3 - diag(Beta)) %*% D
	ans
}

.fs12 = function(X, B,myy, nsamples, pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	ans = (Xsub %*% (( t(Xsub) %*% B)^3)) / dim(Xsub)[2] - 3*B
	ans
}

.fs13 = function(X, B,myy, nsamples, pct.sample)
{
	Ysub = t(X[, .getSamples(nsamples, pct.sample)]) %*% B
	Gpow3 = Ysub^3
	Beta = as.numeric(apply(Ysub %*% Gpow3, 2, "sum"))
	D = diag(1 / (Beta - 3 * dim(t(Ysub))[2]))
	ans = B + myy * B %*% (t(Ysub) %*% Gpow3 - diag(Beta)) %*% D
	ans
}

# tanh
.fs20 = function(X, B, myy, tanh.par, nsamples, pct.sample)
{
	hypTan = tanh(tanh.par * t(X) %*% B)
	ans = X %*% hypTan / nsamples - .ones(dim(B)[1],1) %*% apply(1 - hypTan^2, 2, "sum") * B / nsamples * tanh.par
	ans
}


.fs21 = function(X, B, myy, tanh.par, nsamples, pct.sample)
{
	Y = t(X) %*% B
	hypTan = tanh(tanh.par * Y)
	Beta = apply(Y*hypTan, 2, "sum")
	D = diag(1/(Beta - tanh.par * apply(1 - hypTan^ 2, 2, "sum")))
	ans = B + myy * B %*% (t(Y) %*% hypTan - diag(Beta)) %*% D
	ans
}

.fs22 = function(X, B, myy, tanh.par, nsamples, pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	hypTan = tanh(tanh.par * t(Xsub) %*% B)
	ans = Xsub %*% hypTan / dim(Xsub)[2] - .ones(dim(B)[1],1) %*% apply(1 - hypTan^ 2, 2, "sum") * B / dim(Xsub)[2] * tanh.par
	ans
}

.fs23 = function(X, B, myy, tanh.par, nsamples, pct.sample)
{
	Y = t(X[, .getSamples(nsamples, pct.sample)]) %*% B
	hypTan = tanh(tanh.par * Y)
	Beta = apply(Y * hypTan, 2, "sum")
	D = diag(1 / (Beta - tanh.par * apply(1 - hypTan^ 2, 2, "sum")))
	ans = B + myy * B %*% (t(Y) %*% hypTan - diag(Beta)) %*% D
	ans
}

# gauss
.fs30 = function(X, B, myy, gauss.par, nsamples, pct.sample)
{
	U = t(X) %*% B
	Usquared = U^2
	ex = exp(-gauss.par * Usquared / 2)
	gauss =  U * ex
	dGauss = (1 - gauss.par * Usquared) * ex
	ans = X %*% gauss / nsamples - .ones(dim(B)[1],1) %*% apply(dGauss, 2, "sum") * B / nsamples
	ans
}

.fs31 = function(X, B, myy, gauss.par, nsamples, pct.sample)
{
	Y = t(X) %*% B
	ex = exp(-gauss.par * (Y^2) / 2)
	gauss = Y * ex
	Beta = apply(Y * gauss, 2, "sum")
	D = diag(1 / (Beta - apply( (1 - gauss.par * (Y^2)) * ex, 2, "sum")))
	ans = B + myy * B %*% (t(Y) %*% gauss - diag(Beta)) %*% D
	ans
}


.fs32 = function(X, B, myy, gauss.par, nsamples, pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	U = t(Xsub) %*% B
	Usquared = U^2
	ex = exp(-gauss.par * Usquared / 2)
	gauss =  U * ex
	dGauss = (1 - gauss.par * Usquared) *ex
	ans = Xsub %*% gauss / dim(Xsub)[2] - .ones(dim(B)[1], 1) * apply(dGauss, 2, "sum") * B / dim(Xsub)[2]
	ans
}

.fs33 = function(X, B, myy, gauss.par, nsamples, pct.sample)
{
	Y = t(X[, .getSamples(nsamples, pct.sample)]) %*% B
	ex = exp(-gauss.par * (Y^2) / 2)
	gauss = Y * ex
	Beta = apply(Y * gauss, 2, "sum")
	D = diag(1 / (Beta - sum((1 - gauss.par * (Y^2))* ex)))
	ans = B + myy * B %*% ( t(Y) %*% gauss - diag(Beta)) %*% D				
	ans
}

# skew
.fs40 = function(X, B, myy, nsamples, pct.sample)
{
	ans = (X %*% ((t(X) %*% B)^2)) / nsamples
	ans
}

.fs41 = function(X, B, myy, nsamples, pct.sample)
{
	Y = t(X) %*% B
	Gskew = Y^2
	Beta = apply(Y* Gskew, 2, "sum")
	D = diag(1 / (Beta))
	ans = B + myy * B %*% (t(Y) %*% Gskew - diag(Beta)) %*% D
	ans
}

.fs42 = function(X, B, myy, nsamples, pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	ans = (Xsub %*% ((t(Xsub) %*% B)^2)) / dim(Xsub)[2]
	ans
}

.fs43 = function(X, B, myy, nsamples, pct.sample)
{
	Y = t(X[, .getSamples(nsamples, pct.sample)]) %*% B
	Gskew = Y^2
	Beta = apply(Y* Gskew, 2, "sum")
	D = diag(1 / (Beta))
	ans = B + myy * B %*% (t(Y) %*% Gskew - diag(Beta)) %*% D
	ans
}

#--------------------------------------------------------------------------------------------------------
# deflation methods
#--------------------------------------------------------------------------------------------------------

.fdmethod = function(fmethod, X, w, myy, nsamples, pct.sample, tanh.par, gauss.par)
{
	ans = switch(fmethod,
			f10 = .fd10(X, w,  myy, nsamples,  pct.sample),
			f11 = .fd11(X, w,  myy, nsamples,  pct.sample),
			f12 = .fd12(X, w,  myy, nsamples,  pct.sample),
			f13 = .fd13(X, w,  myy, nsamples,  pct.sample),
			f20 = .fd20(X, w,  myy, tanh.par,  nsamples, pct.sample),
			f21 = .fd21(X, w,  myy, tanh.par,  nsamples, pct.sample),
			f22 = .fd22(X, w,  myy, tanh.par,  nsamples, pct.sample),
			f23 = .fd23(X, w,  myy, tanh.par,  nsamples, pct.sample),
			f30 = .fd30(X, w,  myy, gauss.par, nsamples, pct.sample),
			f31 = .fd31(X, w,  myy, gauss.par, nsamples, pct.sample),
			f32 = .fd32(X, w,  myy, gauss.par, nsamples, pct.sample),
			f33 = .fd33(X, w,  myy, gauss.par, nsamples, pct.sample),
			f40 = .fd40(X, w,  myy, nsamples,  pct.sample),
			f41 = .fd41(X, w,  myy, nsamples,  pct.sample),
			f42 = .fd42(X, w,  myy, nsamples,  pct.sample),
			f43 = .fd43(X, w,  myy, nsamples,  pct.sample))
	return(ans)
}

.fd10 = function(X, w, myy, nsamples, pct.sample)
{
	ans = (X %*% ((t(X) %*% w)^3)) / nsamples - 3 * w
	ans
}

.fd11 = function(X, w, myy, nsamples, pct.sample)
{
	EXGpow3 = (X %*% ((t(X) %*% w)^3)) / nsamples
	Beta = as.numeric(t(w) %*% EXGpow3)
	ans = w - myy * (EXGpow3 - Beta * w) / (3 - Beta)
	ans
}

.fd12 = function(X, w, myy, nsamples, pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	ans = (Xsub %*% ((t(Xsub) %*% w)^3)) / dim(Xsub)[2] - 3 * w
	ans
}

.fd13 = function(X, w, myy, nsamples, pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	EXGpow3 = (Xsub %*% ((t(Xsub) %*% w)^3))/dim(Xsub)[2]
	Beta = as.numeric(t(w)%*%EXGpow3)
	ans = w - myy * (EXGpow3 - Beta * w) / (3 - Beta)
	ans
}

.fd20 = function(X, w,  myy, tanh.par,  nsamples, pct.sample)
{
	hypTan = tanh(tanh.par * t(X) %*% w)
	ans = (X %*% hypTan - tanh.par * sum(1 - hypTan^2) * w) / nsamples
	ans
}

.fd21 = function(X, w,  myy, tanh.par,  nsamples, pct.sample)
{
	hypTan = tanh(tanh.par * t(X) %*% w)
	Beta = as.numeric(t(w) %*% X %*% hypTan)
	ans = w - myy * ((X %*% hypTan - Beta * w)/(tanh.par * sum(1 - hypTan^2) - Beta))
	ans
}

.fd22 = function(X, w,  myy, tanh.par,  nsamples, pct.sample)
{
	Xsub =  X[, .getSamples(nsamples, pct.sample)]
	hypTan = tanh(tanh.par * t(Xsub) %*% w)
	ans = (Xsub %*% hypTan - tanh.par * sum(1 - hypTan^2) * w)/dim(Xsub)[2]
	ans
}

.fd23 = function(X, w,  myy, tanh.par,  nsamples, pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	hypTan = tanh(tanh.par * t(Xsub) %*% w)
	Beta = t(w) %*% Xsub %*% hypTan
	ans = w - myy * ((Xsub %*% hypTan - Beta * w)/(tanh.par * sum(1-hypTan^2) - Beta))
	ans
}

.fd30 = function(X, w,  myy, gauss.par, nsamples, pct.sample)
{
	u = t(X) %*% w
	u2 = u^2
	ex = exp(-gauss.par * u2/2)
	gauss =  u*ex
	dGauss = (1 - gauss.par * u2) *ex
	ans = (X %*% gauss - sum(dGauss) * w) / nsamples
	ans
}

.fd31 = function(X, w,  myy, gauss.par, nsamples, pct.sample)
{
	u = t(X) %*% w
	u2 = u^2
	ex = exp(-gauss.par * u2/2)
	gauss =  u*ex
	dGauss = (1 - gauss.par * u2)*ex
	Beta = as.numeric(t(w) %*% X %*% gauss)
	ans = w - myy * ((X %*% gauss - Beta * w) / (sum(dGauss) - Beta))
	ans
}

.fd32 = function(X, w,  myy, gauss.par, nsamples, pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	u = t(Xsub) %*% w
	u2 = u^2
	ex = exp(-gauss.par * u2/2)
	gauss =  u*ex
	dGauss = (1 - gauss.par * u2) *ex
	ans = (Xsub %*% gauss - sum(dGauss) * w)/dim(Xsub)[2]
	ans
}

.fd33 = function(X, w,  myy, gauss.par, nsamples, pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	u = t(Xsub) %*% w
	u2 = u^2
	ex = exp(-gauss.par * u2/2)
	gauss =  u*ex
	dGauss = (1 - gauss.par * u2) *ex
	Beta = as.numeric(t(w) %*% Xsub %*% gauss)
	ans = w - myy * ((Xsub %*% gauss - Beta * w)/(sum(dGauss) - Beta))
	ans
}

.fd40 = function(X, w,  myy, nsamples,  pct.sample)
{
	ans = (X %*% ((t(X) %*% w)^2)) / nsamples
	ans
}

.fd41 = function(X, w,  myy, nsamples,  pct.sample)
{
	EXGskew = (X %*% ((t(X) %*% w)^2))/nsamples
	Beta = as.numeric(t(w) %*% EXGskew)
	ans = w - myy * (EXGskew - Beta * w)/(-Beta)
	ans
}

.fd42 = function(X, w,  myy, nsamples,  pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	ans = (Xsub %*% ((t(Xsub) %*% w)^2))/dim(Xsub)[2]
	ans
}

.fd43 = function(X, w,  myy, nsamples,  pct.sample)
{
	Xsub = X[, .getSamples(nsamples, pct.sample)]
	EXGskew = (Xsub %*% ((t(Xsub) %*% w)^2))/dim(Xsub)[2]
	Beta = as.numeric(t(w) %*% EXGskew)
	ans = w - myy * (EXGskew - Beta * w)/(-Beta)
	ans
}
################################################################################

# *****************************************************************
# Copyright (c) Erik G. Learned-Miller, 2004.
# *****************************************************************
# RADICAL   Solve the ICA problem in arbitrary dimension.
#
#    Version 1.1. Major bug fix. Faster entropy estimator.
# 
#    Apr.1, 2004. Major bug fix. Whitening matrix was wrong. Thanks
#      to Sergey Astakhov for spotting this one.
#
#    Mar.28, 2004. Sped up inner loop by about 20# with a better
#      entropy estimator.
#
#    Version 1.0. First release.
#
# R translation-Alexios Ghalanos and code

radical = function(X, n.comp = dim(X)[2], demean = TRUE,
		pca.cov = c("ML", "LW", "ROB", "EWMA"), k = 150, 
		augment = FALSE, replications = 30, sd = 0.175, firstEig = 1, 
		lastEig = dim(X)[1], pcaE = NULL, pcaD = NULL, whiteSig = NULL, whiteMat = NULL, 
		dewhiteMat = NULL, rseed = NULL, trace = FALSE, ...)
{
	tic = Sys.time()
	# force evaluation of n.comp before transposing else will lazily evaluate
	# afterwards and we'll get the wrong result
	n.comp = n.comp
	# transpose dataset (as expected by most ICA algorithms)
	X = t(X)
	idx = match(tolower(pca.cov[1]), tolower(c("ML", "LW", "ROB", "EWMA")))
	if(is.na(idx)) stop("\nradical-->error: pca.cov choice not recognized.")
	pca.cov = c("ML", "LW", "ROB", "EWMA")[idx]
	
	########################################
	# STAGE 1: Demean Data
	
	# Remove the mean and check the data
	if(demean){ 
		mixedmean = apply(X, 1, "mean")
		mixedsig  = t(apply(X, 1, FUN = function(x) x - mean(x) ))
	} else{
		mixedmean = rep(0, dim(X)[1])
		mixedsig  = X
	}
	Dim = dim(mixedsig)[1]
	# n.comps + firstEig/lastEig choice
	if(n.comp>nrow(X)) stop("\nCannot Choose more components than signals!")
	if(is.null(firstEig)) firstEig = 1
	if(is.null(lastEig)) lastEig = n.comp
	if(!is.null(firstEig)){
		if(is.null(lastEig)){
			lastEig = firstEig+n.comp
			if(lastEig>nrow(X)) stop("\nfirstEig and n.comp combination greater than no. of signals.")
		} else{
			if(lastEig>nrow(X)) stop("\nlastEig greater than no. of signals.")
			n.comp = lastEig - firstEig + 1
			if(n.comp<1) stop("\nCannot have zero factors in PCA!")
		}
	}
	if(!is.null(lastEig)){
		n.comp = lastEig - firstEig + 1
	}
	# PCA E (m x f)
	if(!is.null(pcaE)){
		if(is.null(pcaD)) stop("\npcaD required when passing pcaE.")
		if(!is.matrix(pcaE)) stop("\npcaE must be an m by f matrix (see documentation)")
		if(nrow(pcaE)!=nrow(X)) stop("\npcaE no.rows must equal no.signals.")
		if(ncol(pcaE)!=n.comp) stop("\npcaE no.cols must equal no.comps.")
	}
	# PCA D (f by f)
	if(!is.null(pcaD)){
		if(is.null(pcaE)) stop("\npcaE required when passing pcaD.")
		if(!is.matrix(pcaD)) stop("\npcaD must be an f by f diagonal matrix (see documentation)")
		if(nrow(pcaD)!=n.comp) stop("\npcaD no.rows must equal no.comps.")
		if(ncol(pcaD)!=n.comp) stop("\npcaD no.cols must equal no.comps.")
		# recover the covariance matrix (for non-dimensionality reduced systems)
		if(NROW(pcaE) == NCOL(pcaE)) C = pcaE %*% pcaD %*% solve(pcaE) else C = NULL
	}
	
	# whiteSig (n by f)
	if(!is.null(whiteSig)){
		if(!is.matrix(whiteSig)) stop("\nwhiteSig must be an n by f matrix (see documentation)")
		if(nrow(whiteSig)!=ncol(X)) stop("\nwhiteSig no.rows must equal no.rows Data.")
		if(ncol(whiteSig)!=n.comp) stop("\nwhiteSig no.cols must equal no.comps.")
		# now transpose
	} else{
		if(!is.null(pcaE)){
			whiteMat = solve(sqrt(pcaD))%*%t(pcaE)
			whiteSig = t(mixedsig) %*% t(whiteMat)
			whiteSig = t(whiteSig)
		}
	}
	
	if(!is.null(dewhiteMat)){
		if(is.null(whiteMat)) stop("\ndewhiteMat cannot be provided without whiteMat.")
		if(nrow(dewhiteMat)!=nrow(X)) stop("\ndewhiteMat must be an m by f matrix (see documentation)")
		if(ncol(dewhiteMat)!=n.comp) stop("\ndewhiteMat must be an m by f matrix (see documentation)")
	}
	
	# whiteMat (m by f)
	if(!is.null(whiteMat)){
		if(!is.matrix(whiteMat)) stop("\nwhiteMat must be an f by m matrix (see documentation)")
		if(nrow(whiteMat)!=n.comp) stop("\nwhiteMat no.rows must equal no.comps.")
		if(ncol(whiteMat)!=nrow(X)) stop("\nwhiteMat no.cols must equal no.signals.")
		if(is.null(dewhiteMat)){
			if(ncol(whiteMat)==nrow(whiteMat)){
				dewhiteMat = t(solve(whiteMat))
			} else{
				stop("\ndewhiteMat needs to be supplied for user supplied non-symmetric whiteMat")
			}
		}
		# now transpose
		if(is.null(whiteSig)){
			whiteSig = t(mixedsig) %*% t(whiteMat)
			whiteSig = t(whiteSig)
		}
	}
	
	if(is.null(rseed)) rseed = runif(1, 1, 100000)
	# When augment is FALSE, do not augment data. Use original data only.
	if(!augment) replications = 1
	Dim = dim(X)[1]
	nsamples = dim(X)[2]
	# m for use in m-spacing estimator.
	m = floor(sqrt(nsamples))
	# Whitening
	# PCA
	jumpPCA = 0
	E = NULL
	D = NULL
	whiteningsig = NULL
	whiteningMatrix = NULL
	dewhiteningMatrix = NULL
	if(!is.null(pcaE))
	{
		jumpPCA = jumpPCA + 1
		E = pcaE
	}
	if(!is.null(pcaD))
	{
		jumpPCA = jumpPCA + 1
		D = pcaD
	}
	# calculate if there are enought parameters to skip PCA and whitening
	jumpWhitening = 0
	if(!is.null(whiteSig))
	{
		jumpWhitening = jumpWhitening + 1
		whiteningsig = whiteSig
	}
	if(!is.null(whiteMat))
	{
		jumpWhitening = jumpWhitening + 1
		whiteningMatrix = whiteMat
	}
	if(!is.null(dewhiteMat))
	{
		jumpWhitening = jumpWhitening + 1
		dewhiteningMatrix = dewhiteMat
	}
	if(trace){
		cat(paste("\nNo. of signals: ", Dim, sep = ""))
		cat(paste("\nNo. of samples: ", nsamples, "\n", sep = ""))
	}
	if(Dim > nsamples && trace)
		warning("\radical-->warning: The signal matrix(X) may be orientated in the wrong way...\n")
	if(jumpWhitening == 3){
		if(trace) cat("\nPCA inputs supplied. PCA calculations not needed.\n")
	} else{
		if(jumpPCA == 2) {
			if(trace) cat("\nValues for PCA calculations supplied. PCA calculations not needed.\n")
		} else{
			if(jumpPCA == 1 && trace) cat("\nYou must supply both pcaE and pcaD in order to jump PCA.\n")
			tmp = .pca(mixedsig, firstEig, lastEig, pca.cov, trace, ...)
			E = tmp$E
			D = tmp$D
			C = tmp$C
		}
	}
	if(jumpWhitening == 3){
		if(trace) cat("\nWhitening not needed.\n")
	} else{
		tmp = whitenv(mixedsig, E, D, trace)
		whiteningsig = tmp$newVectors
		whiteningMatrix = tmp$whiteningMatrix
		dewhiteningMatrix = tmp$dewhiteningMatrix
	}
	Dim = dim(whiteningsig)[1]
	sweeps = Dim - 1
	oldTotalRot = diag(Dim)
	# Current sweep number.
	sweepIter = 0
	totalRot = diag(Dim)
	xcur = whiteningsig
	
	# K represents the number of rotations to examine on the FINAL
	# sweep. To optimize performance, we start with a smaller number of
	# rotations to examine. Then, we increase the
	# number of angles to get better resolution as we get closer to the
	# solution. For the first half of the sweeps, we use a constant
	# number for K. Then we increase it exponentially toward the finish.
	finalK = k
	startKfloat = (finalK/1.3^(ceiling(sweeps/2)))
	newKfloat = startKfloat
	for(sweepNum in 1:sweeps){
		if(trace) cat(paste("\nSweep: ", sweepNum,"/", sweeps, "\n", sep = ""))
		range = pi/2
		# Compute number of angle samples for this sweep.	
		if(sweepNum >(sweeps/2)){
			newKfloat = newKfloat*1.3
			newK = floor(newKfloat)
		} else{
			newKfloat = startKfloat
			newK = max(30,floor(newKfloat))
		}
		for(i in 1:(Dim-1)){
			for(j in (i+1):Dim){
				if(trace) cat(paste("Unmixing dimensions ", i, "...", j, sep = ""))
				# **********************************************
				# Extract dimensions (i,j) from the current data.
				# **********************************************
				curSubSpace = rbind(xcur[i, ], xcur[j, ])
				
				# ***************************************************
				# Find the best angle theta for this Jacobi rotation.
				# ***************************************************
				tmp = .radicalOptThetaC(X = curSubSpace, sd, m, replications, K = newK, range, rseed, trace)
				thetaStar = tmp$thetaStar
				rotStar = tmp$rotStar
				# *****************************************
				# Incorporate Jacobi rotation into solution.
				# *****************************************
				newRotComponent = diag(Dim)
				newRotComponent[i, i] =  cos(thetaStar)
				newRotComponent[i, j] = -sin(thetaStar)
				newRotComponent[j, i] =  sin(thetaStar)
				newRotComponent[j, j] =  cos(thetaStar)
				totalRot = newRotComponent%*%totalRot
				xcur = totalRot%*%whiteningsig
			}
		}
		oldTotalRot = totalRot
	}
	W = totalRot%*%whiteningMatrix
	S = W %*% mixedsig + (W %*% mixedmean) %*% .ones(1, nsamples)
	U = totalRot
	A = t(t(solve(U)) %*% t(dewhiteningMatrix))
	
	toc = Sys.time() - tic
	
	return(list(A = (A), W = (W), U = t(U), S = t(S), E = E, D = D, C = C, mu = mixedmean, 
					whiteningMatrix = (whiteningMatrix), 
					dewhiteningMatrix = (dewhiteningMatrix), 
					rseed = tmp$rseed, elapsed = toc))
}

# *****************************************************************
# Copyright (c) Erik G. Learned-Miller, 2003.
# *****************************************************************
# R/C translation Alexios Ghalanos 2012

# C++ Version
.radicalOptThetaC = function(X, sd, m, replications, K, range, rseed, trace)
{
	# m is the number of intervals in an m-spacing
	# reps is the number of points used in smoothing
	# K is the number of angles theta to check for each Jacobi rotation.
	d = dim(X)[1]
	N = dim(X)[2]
	# This routine assumes that it gets whitened data.
	# First, we augment the points with reps near copies of each point.		
	if(replications == 1){
		xAug = X
	} else{
		set.seed(rseed)
		xAug = matrix(rnorm(d*N*replications, mean = 0, sd = sd), nrow = d, 
				ncol = N*replications, byrow = TRUE) + .repmat(X, 1, replications)
	}
	# Then we rotate this data to various angles, evaluate the sum of 
	# the marginals, and take the min.		
	#perc = range/(pi/2)
	#numberK = perc*K
	#start = floor(K/2-numberK/2)+1
	#endPt = ceiling(K/2+numberK/2)
	ent = .Call("radicalrot", X = as.matrix(xAug), idx = as.integer(c(m,K)), 
			PACKAGE = "rmgarch")
	tmp = sort.int(ent, index.return = TRUE)
	val = tmp$x
	ind = tmp$ix
	thetaStar= (ind[1]-1)/(K-1)*pi/2-pi/4
	if(trace) cat(paste(" rotated ", round(thetaStar/(2*pi)*360, 2), " degrees.\n", sep = ""))
	rotStar = rbind(c(cos(thetaStar), -sin(thetaStar)), c(sin(thetaStar), cos(thetaStar)))
	return(list(thetaStar = thetaStar, rotStar = rotStar))
}

.vasicekm = function(v,m)
{
	len = length(v)
	vals = sort(v)
	# Note that the intervals overlap for this estimator.
	h = sum( log( vals[(m+1):len] - vals[1:(len-m)] ) )
	return(h)
}
###############################################################################

# adapted from octave orth function of Kurt Hornik
.orthogonal = function(x)
{
	tmp = svd(x)
	u = tmp$u
	s = tmp$d
	m = dim(x)[1]
	n = dim(x)[2]
	if(n == 1 || m == 1)  s = s[1] else s = s
	tol = max(m,n) * s[1] * .Machine$double.eps
	rnk = sum(s > tol)
	ans = u[,1:rnk]
	return(ans)
}

# matrix square root via svd
.sqrtm = function(x)
{
	tmp = svd(x)
	sqrtx = tmp$u%*%sqrt(diag(tmp$d))%*%t(tmp$u)
	return(sqrtx)
}

.getSamples = function(max.n, percentage)
{
	ans = numeric()
	while(length(ans) == 0){
		ans = which(runif(max.n) < percentage)
	}
	return(ans)
}

.pca = function(vectors, firstEig = 1, lastEig = dim(vectors)[1], pca.cov, trace, ...)
{
	oldDimension = dim(vectors)[1]
	y = t(vectors)
	covarianceMatrix <- switch(pca.cov,
			ML = (t(y) %*% y)/dim(y)[1],
			LW = lw.cov(y, demean = FALSE,...),
			ROB = cov.rob(y, ...)$cov,
			EWMA = ewma.cov(y, demean = FALSE, ...))
	tmp = eigen(covarianceMatrix)
	D = diag(tmp$values)
	E = tmp$vectors
	rankTolerance = 1e-7
	maxLastEig = sum(diag(D) > rankTolerance)
	if(maxLastEig == 0)
	{
		cat("\nEigenvalues of the covariance matrix are all smaller than tolerance
						of 1e-7. Please make sure that your data matrix contains. nonzero values. If the values 
						are very small try rescaling the data matrix.\n")
		stop("\nafrica-->error: Unable to continue, aborting.\n", call. = FALSE)
	}
	eigenvalues = sort(diag(D), decreasing  = TRUE)
	if(lastEig > maxLastEig)
	{
		lastEig = maxLastEig
		if(trace) cat(paste("\nDimension reduced to ",lastEig-firstEig+1," due to the singularity of covariance matrix\n", sep = ""))
	} else{
		if(trace){
			if(oldDimension == (lastEig - firstEig + 1)){
				cat("Dimension not reduced.\n")
			} else{
				cat("Reducing dimension...\n")
			}
		}
	}
	# Drop the smaller eigenvalues
	if(lastEig < oldDimension)
	{
		lowerLimitValue = (eigenvalues[lastEig] + eigenvalues[lastEig + 1])/ 2
	} else{
		lowerLimitValue = eigenvalues[oldDimension] - 1
	}
	lowerColumns = diag(D) > lowerLimitValue
	# Drop the larger eigenvalues
	if(firstEig > 1){
		higherLimitValue = (eigenvalues[firstEig - 1] + eigenvalues[firstEig]) / 2
	} else{
		higherLimitValue = eigenvalues[1] + 1
	}
	higherColumns = diag(D) < higherLimitValue
	
	# Combine the results from above
	selectedColumns = lowerColumns & higherColumns
	
	# print some info for the user
	if(trace) cat(paste("Selected ", sum (selectedColumns), "dimensions.\n", sep = ""))
	
	if(sum(selectedColumns) != (lastEig - firstEig + 1)) 
		stop("\nafrica-->error: Selected wrong number of dimensions.\n", call. = FALSE)
	
	if(trace)
	{
		cat(paste("Smallest remaining (non-zero) eigenvalue ", round(eigenvalues[lastEig],5), "\n", sep = ""))
		cat(paste("Largest remaining (non-zero) eigenvalue ", round(eigenvalues[firstEig],5), "\n", sep = ""))
		cat(paste("Sum of removed eigenvalues ", round(sum(diag(D)*(!selectedColumns)), 5), "\n", sep = ""))
	}
	
	# Select the colums which correspond to the desired range
	# of eigenvalues.
	E = .selcol(E, selectedColumns)
	D = .selcol(t(.selcol(D, selectedColumns)), selectedColumns)
	if(trace)
	{
		sumAll = sum(eigenvalues)
		sumUsed = sum(diag(D))
		retained = (sumUsed / sumAll) * 100
		cat(paste(round(retained,2), " % of (non-zero) eigenvalues retained.\n", sep = ""))
	}
	return(list(E = E, D = D, C = covarianceMatrix))
}


.selcol = function(oldMatrix, maskVector)
{
# newMatrix = selcol(oldMatrix, maskVector);
# Selects the columns of the matrix that marked by one in the given vector.
# The maskVector is a column vector.
	takingMask = numeric()
	if(length(maskVector) != dim(oldMatrix)[2])
		stop("\nThe mask vector and matrix are of uncompatible size.\n", call. = FALSE)
	numTaken = 0
	for(i in 1:length(maskVector))
	{
		if(maskVector[i] == 1)
		{
			takingMask[numTaken + 1] = i
			numTaken = numTaken + 1
		}
	}
	newMatrix = oldMatrix[, takingMask]
	return(newMatrix)
}

# A translation of the matlab function from the FASTICA toolbox
# http://www.cis.hut.fi/projects/ica/fastica/

whitenv  = function(vectors, E, D, trace)
{
	if (any(diag(D) < 0))
		stop("\nafrica-->error: negative eigenvalues computed from the 
						covariance matrix. These are due to rounding errors 
						in R (the correct eigenvalues are probably very small). 
						To correct the situation please reduce the number of 
						dimensions in the data by using the lastEig argument in 
						function the ICA function.\n", call. = FALSE)
	whiteningMatrix = solve(sqrt(D))%*%t(E)
	dewhiteningMatrix = E %*% sqrt(D)
	if(trace) cat("Whitening...\n")
	# newVectors calculated for check
	newVectors =  whiteningMatrix %*% vectors
	if(any(is.complex(newVectors))) stop("\nfastica-->error: Whitened vectors have imaginary values.\n", call. = FALSE)
	
	return(list(newVectors = newVectors, whiteningMatrix = whiteningMatrix, 
					dewhiteningMatrix = dewhiteningMatrix))
}
.zeros = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(0, nrow = n, ncol = m)
	return(sol)
}

.ones = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(1, nrow = n, ncol = m)
	return(sol)
}

# frobenius norm
.fnorm = function(x)
{
	sqrt(sum(as.numeric(x)^2))
}

# matrix rotation
.rot90  = function(a, n=1)
{
	n <- n %% 4
	if (n > 0) {
		return (Recall( t(a)[nrow(a):1,],n-1) )
	} else {
		return (a)
	}
}

.repmat = function(a,n,m) {kronecker(matrix(1,n,m),a)}



# Ledoit and Wolfe Shrinkage Estimator
lw.cov = function(X, shrink = -1, demean = FALSE){
	n = dim(X)[1]
	m = dim(X)[2]
	meanx = colMeans(X)
	if(demean) X = scale(X, scale = FALSE)
	# compute sample covariance matrix
	samplecov = (t(X)%*%X)/n
	# compute prior
	meanvar = mean(diag(samplecov))
	prior = meanvar * diag(m)
	if(shrink == -1){
		# compute shrinkage parameters
		# p in paper 
		y = X^2
		phiMat = (t(y)%*%y)/n - 2*(t(X)%*%X)*samplecov/n + samplecov^2
		phi = sum(apply(phiMat, 1, "sum"))
		# c in paper
		cgamma = norm(samplecov-prior,'F')^2
		# shrinkage constant
		kappa = phi/cgamma
		shrinkage = max(0, min(1, kappa/n))
	} else{
		shrinkage = shrink
	}
	sigma = shrinkage * prior + (1 - shrinkage)*samplecov
	return(sigma)
}

ewma.cov = function(X, lambda = 0.96, demean = FALSE)
{
	n = dim(X)[1]
	i = (0:(n-1))
	ewma.wt = lambda^i
	ewma.wt = ewma.wt/sum(ewma.wt)
	covm = cov.wt(X, wt = rev(ewma.wt), center = demean)$cov
	return(covm)
}

# check + reconstruct demo:
# Data = dji30ret[,1:4]
check_ICA<-function(Data, A, W, U, E, D, S, mu, whiteningMatrix, dewhiteningMatrix){
	cat("\nICA system checks...")
	cat("\nOriginal Data:")
	cat(paste("\nn.rows (n): ", dim(Data)[1], sep = ""))
	cat(paste("\nn.cols (m): ", dim(Data)[2], sep = ""))
	cat(paste("\nn.comp (f): ", dim(U)[2], sep = ""))
	cat("\n")
	print(head(Data, 4))
	cat("\nDemean Data (dX) = Data - mu:")
	dX = t(apply(Data, 1, FUN = function(x) x - mu))
	cat("\n")
	print(head(dX, 4))
	cat("\nWhiten Data (PCA):")
	cat("\nwhiteningMatrix [f by m] == solve(sqrt(D))%*%t(E)")
	cat("\n")
	# 1/sqrt(D)
	print(whiteningMatrix == solve(sqrt(D))%*%t(E))
	cat("\nWhitenedData [n by f] = dX %*% t(whiteningMatrix)")
	cat("\n")
	print(head(dX %*% t(whiteningMatrix), 4))
	cat("\nW [f by m] = solve(dewhiteningMatrix %*% U)")
	cat("\nW [f by m] = t(t(whiteningMatrix) %*% U) [for dimensionality reduced systems]\n")	
	cat("\n")
	print(round(W,8) == round(t(t(whiteningMatrix) %*% U),8))
	if(dim(U)[2] == ncol(Data)) print(round(W,8) == round(solve(dewhiteningMatrix %*% U),8))
	cat("\nMixing Matrix (A) [f by m]  == dewhiteningMatrix %*% U")
	cat("\n")
	print(round(A,8)  == round(dewhiteningMatrix %*% U, 8))
	cat("\nUnconditional Covariance Matrix  == A %*% t(A)")
	cat("\nt(dX)%*%dX/n == A %*% t(A) [Not true for dimensionality reduced system]\n")
	print(round(t(dX)%*%dX/NROW(dX),8) == round(A %*% t(A),8))
	cat("\n")
	cat("\nIndependent Factors S [n by f] == dX %*% W + matrix(mu %*% W, ncol = dim(dX)[2], nrow = dim(dX)[1], byrow = TRUE)")
	cat("\n")
	Fc = dX %*% t(W)
	print(head(Fc + matrix(mu %*% t(W), ncol = dim(Fc)[2], nrow = dim(Fc)[1], byrow = TRUE), 5) == head(S, 5))
	cat("\nOriginal Data (Data) == S %*% t(A) [Not true for dimensionality reduced system]")
	cat("\n")
	print(round(head(S %*% t(A), 4),8) == round(head(Data, 4),8))
	cat("\nend checks\n")
	return(invisible(0))
}
##################################################################################
