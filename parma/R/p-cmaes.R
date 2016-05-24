#################################################################################
##
##   R package parma
##   Alexios Ghalanos Copyright (C) 2012-2013 (<=Aug)
##   Alexios Ghalanos and Bernhard Pfaff Copyright (C) 2013- (>Aug)
##   This file is part of the R package parma.
##
##   The R package parma is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package parma is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

################################################################################
# This file contains my translation into R of Hansen's cmaes solver version 3.6
# http://www.lri.fr/~hansen/cmaes_inmatlab.html#matlab

# Currently, the following functionality is NOT implemented
# 1. noisy function evalution model
# 2. resume from file functionality
# 3. plotting
################################################################################

min.int = function(x, index.return = FALSE){
	idx = which(x == min(x))[1]
	x = x[idx]
	if(index.return){
		return(list(x = x, ix = idx))
	} else{
		return(x)
	}
}

max.int = function(x, index.return = FALSE){
	idx = which(x == max(x))[1]
	x = x[idx]
	if(index.return){
		return(list(x = x, ix = idx))
	} else{
		return(x)
	}
}
# some translation functions taken from the pracma package

find <- function(v){ which( if (is.logical(v)) v else v != 0 ) }

sqrtm <- function(A, kmax = 20, tol = .Machine$double.eps^(1/2)) {
	stopifnot(is.numeric(A), is.matrix(A))
	if (nrow(A) != ncol(A))
		stop("Matrix 'A' must be square.")
	
	# should be "try(solve(A))"
	P0 <- A; Q0 <- diag(nrow(A))
	
	k <- 0  # then k <- 1
	while (norm(A - P0 %*% P0, 'F') > tol && k < kmax) {
		P1 <- 0.5 * (P0 + solve(Q0))
		Q1 <- 0.5 * (Q0 + solve(P0))
		P0 <- P1
		Q0 <- Q1
		k <- k + 1
	}
	
	# k < 0 if iteration has not converged.
	if (k >= kmax) k <- -1
	
	# return sqrtm(A) and sqrtm(A)^-1
	return(list(B = P0, Binv = Q0, k = k, acc = norm(A - P0 %*% P0, 'F')))
}



eye <- function(n, m = n) {
	stopifnot(is.numeric(n), length(n) == 1,
			is.numeric(m), length(m) == 1)
	n <- floor(n)
	m <- floor(m)
	if (n <= 0 || m <= 0) return(matrix(NA, 0, 0))
	else                  return(base::diag(1, n, m))
}

ones <- function(n, m = n) {
	stopifnot(is.numeric(n), length(n) == 1,
			is.numeric(m), length(m) == 1)
	n <- floor(n)
	m <- floor(m)
	if (n <= 0 || m <= 0) return(matrix(1, 0, 0))
	else                  return(matrix(1, n, m))
}

zeros <- function(n, m = n) {
	stopifnot(is.numeric(n), length(n) == 1,
			is.numeric(m), length(m) == 1)
	n <- floor(n)
	m <- floor(m)
	if (n <= 0 || m <= 0) return(matrix(0, 0, 0))
	else                  return(matrix(0, n, m))
}

tril <- function(M, k = 0) {
	if (k == 0) {
		M[upper.tri(M, diag = FALSE)] <- 0
	} else {
		M[col(M) >= row(M) + k + 1] <- 0
	}
	return(M)
}


triu <- function(M, k = 0) {
	if (k == 0) {
		M[lower.tri(M, diag = FALSE)] <- 0
	} else {
		M[col(M) <= row(M) + k - 1] <- 0
	}
	return(M)
}


repmat <- function(a, n, m = n) {
	if (length(a) == 0) return(c())
	if (!is.numeric(a) && !is.complex(a))
		stop("Argument 'a' must be a numeric or complex.")
	if (is.vector(a))
		a <- matrix(a, nrow = 1, ncol = length(a))
	if (!is.numeric(n) || !is.numeric(m) ||
			length(n) != 1 || length(m) != 1)
		stop("Arguments 'n' and 'm' must be single integers.")
	n <- max(floor(n), 0)
	m <- max(floor(m), 0)
	if (n <= 0 || m <= 0)
		return(matrix(0, nrow = n, ncol = m))
	
	matrix(1, n, m) %x% a  # Kronecker product
}

reshape <- function(a, n, m) {
	if (missing(m)) m <- length(a) %/% n
	if (length(a) != n*m)
		stop("Matrix 'a' does not have n*m elements")
	dim(a) <- c(n, m)
	return(a)
}

cond <- function(M, p = 2) {
	if (length(M) == 0)
		return(0)
	if (!is.numeric(M))
		stop("Argument 'M' must be a numeric matrix.")
	if (is.vector(M))
		M <- matrix(c(M), nrow = length(M), ncol = 1)
	
	if (length(M) == 0) return(c())
	if (ncol(M) != nrow(M) && p != 2)
		stop("Matrix 'M' must be square if p = 2.")
	
	if (p == 2) {
		s <- svd(M)$d
		cnd <- if (any(s == 0)) Inf else max(s) / min(s)
	} else {
		stop("At the moment, p-norms other than p = 2 are not implemented.")
		#cnd <- norm(M, p) * norm(inv(M), p)
	}
	return(cnd)
}

.eval.control = function(value, control){
	i = which(value == names(control))
	if(length(i)>0){
		tmp = eval.parent(parse(text=paste(control[[i]],sep="")))
	} else{
		tmp = NA
	}
	return(tmp)
}

cmaes.control = function(
		options = list(StopFitness = -Inf, MaxFunEvals = Inf,
				MaxIter = '1e3*(N+5)^2/sqrt(popsize)', StopFunEvals = Inf,
				StopIter = Inf, TolX = '1e-11*max(insigma)', TolUpX = '1e3*max(insigma)',
				TolFun = 1e-12, TolHistFun = 1e-13, StopOnStagnation = TRUE,
				StopOnWarnings = TRUE, StopOnEqualFunctionValues = '2 + N/3',
				DiffMaxChange = Inf, DiffMinChange = 0, WarnOnEqualFunctionValues = FALSE,
				EvalParallel = FALSE, EvalInitialX = TRUE, Restarts = 0,
				IncPopSize = 2, PopSize = '4 + floor(3*log(N))', ParentNumber = 'floor(popsize/2)',
				RecombinationWeights = c("superlinear", "linear", "constant"), 
				DiagonalOnly = '0*(1+100*N/sqrt(popsize))+(N>=1000)', 
				CMA = TRUE, Seed = 'as.integer(Sys.time())', DispFinal = TRUE, 
				DispModulo = 100, Warnings = FALSE), 
		CMA = list(cs = '(mueff+2)/(N+mueff+3)', 
				damps = '1 + 2*max(0,sqrt((mueff-1)/(N+1))-1) + cs',
				ccum = '(4 + mueff/N) / (N+4 + 2*mueff/N)',
				ccov1 = '2 / ((N+1.3)^2+mueff)',
				ccovmu = '2 * (mueff-2+1/mueff) / ((N+2)^2+mueff)',
				active = 0))
{

	xoptions = list(StopFitness = -Inf, MaxFunEvals = Inf,
			MaxIter = '1e3*(N+5)^2/sqrt(popsize)', StopFunEvals = Inf,
			StopIter = Inf, TolX = '1e-11*max(insigma)', TolUpX = '1e3*max(insigma)',
			TolFun = 1e-12, TolHistFun = 1e-13, StopOnStagnation = TRUE,
			StopOnWarnings = TRUE, StopOnEqualFunctionValues = '2 + N/3',
			DiffMaxChange = Inf, DiffMinChange = 0, WarnOnEqualFunctionValues = FALSE,
			EvalParallel = FALSE, EvalInitialX = TRUE, Restarts = 0,
			IncPopSize = 2, PopSize = '4 + floor(3*log(N))', ParentNumber = 'floor(popsize/2)',
			RecombinationWeights = c("superlinear"), 
			DiagonalOnly = '0*(1+100*N/sqrt(popsize))+(N>=1000)', 
			CMA = TRUE, Seed = 'as.integer(Sys.time())', DispFinal = TRUE, 
			DispModulo = 100, Warnings = FALSE)
	xCMA = list(cs = '(mueff+2)/(N+mueff+3)', 
			damps = '1 + 2*max(0,sqrt((mueff-1)/(N+1))-1) + cs',
			ccum = '4/(N+4)',
			ccov1 = '2 / ((N+1.3)^2+mueff)',
			ccovmu = '2 * (mueff-2+1/mueff) / ((N+2)^2+mueff)',
			active = 0)
	
	if(is.null(options)) options = xoptions
	if(is.null(CMA)) CMA = xCMA
	
	idx1 = sapply(names(options), FUN = function(x) match(x, names(unlist(xoptions))))
	if(length(idx1)>0) for(i in 1:length(idx1)) xoptions[[idx1[i]]] = options[[i]][1]
	idx3 = sapply(names(CMA), FUN = function(x) match(x, names(unlist(xCMA))))
	if(length(idx3)>0) for(i in 1:length(idx3)) xCMA[[idx3[i]]] = CMA[[i]][1]
	return(list(options = xoptions, cma = xCMA))
}

myprctile = function(inar, perc, idx = NULL){
	# Computes the percentiles in vector perc from vector inar
	# returns vector with length(res)==length(perc)
	# idx: optional index-array indicating sorted order
	#
	N = length(inar)
	flgtranspose = 0
	#if(dim(perc)[1] > 1){
	#	perc = t(perc)
	#	flgtranspose = 1
	#	if(dim(perc)[1] > 1) stop('\nperc must not be a matrix')
	#
	#}
	#if(dim(inar)[1] > 1 && dim(inar)[2] > 1) stop('\ndata inar must not be a matrix')
	# sort inar
	if(is.null(idx)){
		tmp = sort.int(inar, index.return = TRUE)
		sar = tmp$x
		idx = tmp$ix
	} else{
		sar = inar[idx]
	}
	res = NULL
	for(p in perc){
		if(p <= 100*(0.5/N)){
			res = c(res, sar[1])
		} else if(p >= 100*((N-0.5)/N)){
			res = c(res, sar[N])
		} else{
			# find largest index smaller than required percentile
			availablepercentiles = 100*((1:N)-0.5)/N
			i = max(which(p > availablepercentiles))
			# interpolate linearly
			res = c(res, sar[i]	+(sar[i+1] - sar[i])*(p - availablepercentiles[i])/
							(availablepercentiles[i+1] - availablepercentiles[i]))
		}
	}
	return(res)
}
local_noisemeasurement = function(arf1, arf2, lamreev, theta, cutlimit = Inf){
	arf1dim2 = dim(arf1)[2]
	arf2dim2 = dim(arf2)[2]
	arf1dim1 = dim(arf1)[1]
	arf2dim1 = dim(arf2)[1]
	
	if(arf1dim1 != 1){
		stop('\narf1 must be an 1xlambda matrix')
	}
	if(arf2dim1 != 1){
		stop('\narf1 must be an 1xsomething matrix')
	}
	if(arf1dim2 < arf2dim2){
		# not really necessary, but saver
		stop('\narf2 must not be smaller than arf1 in columns')
	}
	lam = arf1dim2
	if(arf1dim2 != arf2dim2){
		arf2 = cbind(arf2, arf1[,(arf2dim2+1):lam, drop = FALSE])
	}
	# capture unusual values
	diffarf1 = diff(arf1[1,])
	if(any(diffarf1 == 0)){
		# this will presumably interpreted as rank change, because
		# sort(ones(...)) returns 1,2,3,...
		warning(paste("\n",(sum(diffarf1 == 0)), " equal function values", sep=""))
	}
	# compute rank changes into rankDelta
	# compute ranks
	idx = sort.int(c(arf1[1,],arf2[1,]), index.return = TRUE)$ix
	ranks = sort.int(idx, index.return = TRUE)$ix
	ranks = t(matrix(ranks, nrow = lam, ncol = 2))
	rankDelta = matrix(ranks[1, ] - ranks[2,] - sign(ranks[1,] - ranks[2,]), nrow = 1)
	# compute rank change limits using both ranks(1,...) and ranks(2,...)
	sumlim = rep(0, lamreev)
	for(i in 1:lamreev){
		sumlim[i] = 0.5 * (myprctile(abs((1:(2*lam-1) - (ranks[1,i] - (ranks[1,i]>ranks[2,i])))), theta*50) +   
					myprctile(abs((1:(2*lam-1) - (ranks[2,i] - (ranks[2,i]>ranks[1,i])))), theta*50))
	}
	# compute measurement
	#s = abs(rankDelta(1:lamreev)) - sumlim; % lives roughly in 0..2*lambda
	#                               max: 1 rankchange in 2*lambda is always fine
	s = abs(rankDelta[1,(1:lamreev)]) - pmax(1, sumlim) # lives roughly in 0..2*lambda
	# cut-off limit
	idx = abs(s) > cutlimit
	if(length(idx)>0) s[idx] = sign(s[idx]) * cutlimit
	s = mean(s)
	return(list(s = s, ranks = ranks, rankDelta = rankDelta))
}

myrange = function(x){ return( max(x, na.rm = TRUE) - min(x, na.rm = TRUE) ) }

xintobounds = function(x, lbounds, ubounds){
	#
	# x can be a column vector or a matrix consisting of column vectors
	#
	if(is.matrix(x)){
		dimx2 = dim(x)[2]
		dimx1 = dim(x)[1]
	} else{
		dimx2 = 1
		dimx1 = length(x)
	}
	if(!is.null(lbounds)){
		if(length(lbounds)==1){
			idx = x < lbounds
			x[idx] = lbounds
		} else{			
			arbounds = matrix(lbounds, dimx1, dimx2)
		}
		idx = x < arbounds
		x[idx] = arbounds[idx]
	} else{
		idx = 0
	}
	if(!is.null(ubounds)){
		if(length(ubounds)==1){
			idx2 = x > ubounds
			x[idx2] = ubounds
		} else{
			arbounds = matrix(ubounds, dimx1, dimx2)
		}
		idx2 = x > arbounds
		x[idx2] = arbounds[idx2]
	} else{
		idx2 = 0
	}
	idx = idx2 - idx
	return(list(x = x, idx = idx) )
}


cmaes = function(pars, fun, lower = rep(0, length(pars)), upper = rep(1, length(pars)), insigma = 1, ctrl = cmaes.control(), ...){
	# some tests on the validity of fun
	irun=0
	input = list()
	fitness = list()
	bnd = list()
	neg = list()
	timers = list()
	out = list()
	timer = list()
	# some checks on sigma and upper/lower bounds
	counteval = 0
	countevalNA = 0
	opts = ctrl$options
	cma = ctrl$cma
	flg_future_setting = FALSE
	noiseHandling = FALSE
	flgscience = FALSE
	input$fun = fun
	input$pars = pars
	input$sigma = insigma
	Warnings = .eval.control('Warnings', opts)
	if(is.null(insigma)){
		if(all(dim(pars) > 1)){
			insigma = apply(pars,1, sigma)
			if(any(insigma == 0)) stop('\nInitial search volume is zero, choose signa or pars appropriately')
		} else{
			insigma = abs(pars)*2
		}
	}
	arxvalid = arx = arz = NULL
	#insigma = .eval.control('Resume', opts)
	while(irun <= .eval.control('Restarts', opts)){
		cat("\n")
		# for-loop does not work with resume
		irun = irun + 1
		# ------------------------ Initialization -------------------------------
		# Handle resuming of old run
		#flgresume = .eval.control(Resume, opts)
		xmean = pars
		# 
		if(all(dim(as.matrix(xmean)) > 1)){
			xmean = colMeans(xmean)
			N = dim(xmean)[1]
			# in case if xstart is a population
		} else if(!is.null(dim(as.matrix(xmean))) && dim(as.matrix(xmean))[2] > 1){
			xmean = t(xmean)
			N = length(xmean)
		} else{
			N = length(xmean)
		}
		# ignore case of passing pars as matrix of population
		# Assign settings from input parameters and options for myeval...
		numberofvariables = N
		lambda0 = floor(.eval.control('PopSize', opts) * .eval.control('IncPopSize', opts)^(irun-1)) 
		# lambda0 = floor(myeval(opts.PopSize) * 3^floor((irun-1)/2)); 
		popsize = lambda0
		lambda = lambda0
		
		if(!is.null(dim(insigma)) && all(dim(insigma) == c(N,2))){
			insigma = 0.5 * (insigma[,2] - insigma[,1])
		}
		
		stopFitness = .eval.control('StopFitness', opts)
		stopMaxFunEvals = .eval.control('MaxFunEvals', opts)
		stopMaxIter = .eval.control('MaxIter', opts)
		stopFunEvals = .eval.control('StopFunEvals', opts)
		stopIter = .eval.control('StopIter', opts)
		stopTolX = .eval.control('TolX', opts)
		stopTolUpX = .eval.control('TolUpX', opts)
		stopTolFun = .eval.control('TolFun', opts)
		stopTolHistFun = .eval.control('TolHistFun', opts)
		stopOnStagnation = .eval.control('StopOnStagnation', opts) 
		stopOnWarnings = .eval.control('StopOnWarnings', opts) 
		flgreadsignals = .eval.control('ReadSignals', opts)
		flgWarnOnEqualFunctionValues = .eval.control('WarnOnEqualFunctionValues', opts)
		flgEvalParallel = .eval.control('EvalParallel', opts)
		stopOnEqualFunctionValues = .eval.control('StopOnEqualFunctionValues', opts)
		arrEqualFunvals = matrix(0, 1, 10+N)
		flgDiagonalOnly = .eval.control('DiagonalOnly', opts) 
		flgActiveCMA = .eval.control('CMA', opts) 
		#noiseHandling = .eval.control('Noise', opts)
		#noiseMinMaxEvals = .eval.control('minmaxevals', Noise)
		#noiseAlphaEvals = .eval.control('alphaevals', Noise)
		#noiseCallback = .eval.control('callback', Noise) 
		flgdisplay = .eval.control('DispFinal', opts)
		#flgplotting = .eval.control('LogPlot', Noise)
		verbosemodulo = .eval.control('DispModulo',opts)
		# for now avoid the resum
		flgsaving = 0
		flgsavingfinal = 0
		maxdx = .eval.control('DiffMaxChange', opts)
		# maximal sensible variable change
		mindx = .eval.control('DiffMinChange', opts)
		# minimal sensible variable change 
		# can both also be defined as Nx1 vectors
		lbounds = lower
		ubounds = upper
		if(length(lbounds) == 1){
			lbounds = rep(lbounds, N)
		}
		if(length(ubounds) == 1){
			ubounds = rep(ubounds, N)
		}
		if(is.null(insigma)){
			# last chance to set insigma
			if(all(lbounds > -Inf) && all(ubounds < Inf) ){
				if(any(lbounds>=ubounds)) stop("\nupper bound must be greater than lower bound")
				insigma = 0.3*(ubounds-lbounds)
				stopTolX = .eval.control('TolX', opts)
				# reevaluate these
				stopTolUpX = .eval.control('TolUpX', opts)
			} else{
				stop("\nInitial step sizes (sigma) not determined")
			}
		}
		# Check all vector sizes
		if(dim(as.matrix(xmean))[2] > 1 || dim(as.matrix(xmean))[1] != N){
			stop(paste("\nintial search point should be a vector of size ", N, sep = ""))
		} 
		if(length(insigma) !=1 && length(insigma)!=N){
			stop(paste("\ninput parameter sigma should be a scalar or a vector of size ", N, sep = ""))
		} 
		if(length(stopTolX)>1 && length(stopTolX)!=N){
			stop(paste("\noption TolX should be a scalar or a vector of size ", N, sep=""))
		} 
		if(length(stopTolUpX)>1 && length(stopTolUpX)!=N){
			stop(paste("\noption TolUpX should be a scalar or a vector of size ", N, sep=""))
		} 
		if(length(maxdx)>1 && length(maxdx)!=N){
			stop(paste("\noption DiffMaxChange should be a scalar or a vector of size ", N, sep=""))
		} 
		if(length(mindx)>1 && length(mindx)!=N){
			stop(paste("\noption DiffMinChange should be a scalar or a vector of size ", N, sep=""))
		} 
		if(length(lbounds)>1 && length(lbounds)!=N){
			stop(paste("\noption lower should be a scalar or a vector of size ", N, sep=""))
		} 
		if(length(ubounds)>1 && length(ubounds)!=N){
			stop(paste("\noption upper should be a scalar or a vector of size ", N, sep=""))
		}
		# Initialize dynamic internal state parameters
		if(any(insigma <= 0)) stop("\nInitial search volume (sigma) must be greater than zero")
		if(max(insigma)/min(insigma) > 1e6) stop("\nInitial search volume (sigma) badly conditioned")
		sigma = max(insigma)
		# overall standard deviation
		pc = ps = rep(0, N)
		# evolution paths for C and sigma
		if(length(insigma) == 1) insigma = rep(insigma, N)
		diagD = insigma/max(insigma)
		# diagonal matrix D defines the scaling
		diagC = diagD^2 
		if(flgDiagonalOnly != 1){
			# use at some point full covariance matrix
			B = diag(N)
			# B defines the coordinate system
			BD = B * repmat(t(diagD), N, 1)
			# B*D for speed up only
			C = diag(diagC)
			# covariance matrix == BD*(BD)'
		}
		if(flgDiagonalOnly == 1){
			B = 1
		}
		fitness$hist = rep(NA, 10+ceiling(3*10*N/lambda))
		# history of fitness values
		fitness$histsel = rep(NA, ncol = 10+ceiling(3*10*N/lambda)) 
		# history of fitness values
		fitness$histbest = NULL
		# history of fitness values
		fitness$histmedian = NULL
		# history of fitness values
		# Initialize boundary handling
		bnd$isactive = (any(lbounds > -Inf) || any(ubounds < Inf) )
		if(bnd$isactive){
			if(any(lbounds>ubounds)){
				stop("\nlower bound found to be greater than upper bound")
			}
			tmp = xintobounds(xmean, lbounds, ubounds)
			xmean = tmp$x
			ti = tmp$idx
			# just in case
			if(length(ti)>0 && any(ti)!=0 && Warnings) warning("\nInitial point was out of bounds, corrected")
			bnd$weights = rep(0, N)
			# weights for bound penalty
			# scaling is better in axis-parallel case, worse in rotated
			bnd$flgscale = 0
			# scaling will be omitted if zero 
			if(bnd$flgscale != 0){ 
				bnd$scale = diagC/mean(diagC)
			} else{
				bnd$scale = rep(1, N)
			}
			idx = as.numeric(lbounds > -Inf | ubounds < Inf)
			if(length(idx) == 1){
				idx = rep(idx, N)
			}
			bnd$isbounded = rep(0, N)
			bnd$isbounded[which(idx != 0)] = 1
			maxdx = pmin(maxdx, (ubounds - lbounds)/2)
			if(any(sigma*sqrt(diagC) > maxdx)){
				fac = min(maxdx/sqrt(diagC))/sigma
				sigma = min(maxdx/sqrt(diagC))
				if(Warnings) warning(paste("Initial sigma multiplied by the factor ",  fac, ", because it was larger than half of one of the boundary intervals", sep = ""))
			}
			idx = as.numeric(lbounds > -Inf & ubounds < Inf)
			dd = diagC
			if(any(5*sigma*sqrt(dd[idx]) < (ubounds[idx] - lbounds[idx]))){
				if(Warnings) warning("\nInitial sigma is, in at least one coordinate	much smaller than the given boundary intervals. For reasonable global search performance sigma\n 
should be between 0.2 and 0.5 of the bounded interval in each coordinate. If all coordinates have lower and upper bounds sigma can be empty.")
			}
			bnd$dfithist = 1
			#delta fit for setting weights
			bnd$aridxpoints = NULL
			# remember complete outside points
			bnd$arfitness = NULL
			# and their fitness
			bnd$validfitval = 0
			bnd$iniphase = 1
		}
		# ooo initial feval, for output only
		if(irun == 1){
			out$solutions$bestever$x = xmean
			out$solutions$bestever$f = Inf
			# for simpler comparison below
			out$solutions$bestever$evals = counteval
			bestever = out$solutions$bestever
		}
		if(.eval.control('EvalInitialX', opts)){
			fitness$histsel[1] = fitness$hist[1] = fun(xmean, ...)
			#fitness$histsel[1] = fitness$hist[1] = fun(xmean, Data, forecast) 
			counteval = counteval + 1
			if(fitness$hist[1] < out$solutions$bestever$f){
				out$solutions$bestever$x = xmean
				out$solutions$bestever$f = fitness$hist[1]
				out$solutions$bestever$evals = counteval
				bestever = out$solutions$bestever
			}
		} else{
			fitness$histsel[1] = fitness$hist[1] = NA
		}
		# initialize random number generator
		rseed = .eval.control('Seed', opts)
		set.seed(rseed)
		#qqq
		#  load(opts.SaveFilename, 'startseed');
		#  randn('state', startseed);
		#  disp(['SEED RELOADED FROM ' opts.SaveFilename]);
		startseed = opts$Seed
		# for retrieving in saved variables
		# Initialize further constants
		chiN = N^0.5*(1-1/(4*N)+1/(21*N^2))
		# expectation of ||N(0,I)|| == norm(randn(N,1))
		countiter = 0
		# Initialize records and output
		if(irun == 1){
			timer$t0 = Sys.time()
			# TODO: keep also median solution? 
			out$evals = counteval
			# should be first entry
			out$stopflag = NULL
			outiter = 0
			# Write headers to output data files
			filenameprefix = opts$LogFilenamePrefix
		}
	######################################################	
		stopflag = NULL
		lambda_hist = NULL
		lambda_last = NA
		while(is.null(stopflag)){
			# set internal parameters
			if(countiter == 0 || lambda != lambda_last){
				if(countiter > 0 && floor(log10(lambda)) != floor(log10(lambda_last)) && flgdisplay){
					cat(paste("\n  lambda = ", lambda, sep = ""))
					lambda_hist = cbind(lambda_hist,  rbind(countiter+1, lambda))
				} else{
					lambda_hist = rbind(countiter+1, as.numeric(lambda)) 
				}
				lambda_last = lambda
				# Strategy internal parameter setting: Selection  
				mu = .eval.control('ParentNumber', opts)
				# number of parents/points for recombination
				if( opts$RecombinationWeights == "equal"){
					weights = rep(1, mu)
					#(mu_I,lambda)-CMA-ES
				} else if(opts$RecombinationWeights == "linear"){
					weights = mu + 0.5 - (1:mu)
				} else{
					# superlinear
					weights = log(mu + 0.5)-log(1:mu)
					# muXone array for weighted recombination
					# qqq mu can be non-integer and
					# should become ceil(mu-0.5) (minor correction)	
				}
				mueff = sum(weights)^2/sum(weights^2)
				# variance-effective size of mu
				weights = weights/sum(weights)
				# normalize recombination weights array
				if(mueff == lambda){
					stop("\nCombination of values for PopSize, ParentNumber and and RecombinationWeights is not reasonable")
				}	
				# Strategy internal parameter setting: Adaptation
				cc = .eval.control('ccum', cma)
				# time constant for cumulation for covariance matrix
				cs = .eval.control('cs', cma)
				# old way TODO: remove this at some point
				# mucov = mueff;   % size of mu used for calculating learning rate ccov
				# ccov = (1/mucov) * 2/(N+1.41)^2 ... % learning rate for covariance matrix
				#        + (1-1/mucov) * min(1,(2*mucov-1)/((N+2)^2+mucov)); 
				# new way
				if(.eval.control('CMA', opts)){
					ccov1 = .eval.control('ccov1', cma)
					ccovmu = min(1 - ccov1, .eval.control('ccovmu', cma))
				} else{
					ccov1 = 0
					ccovmu = 0
				}
				# flgDiagonalOnly = -lambda*4*1/ccov; % for ccov==1 it is not needed
				# 0 : C will never be diagonal anymore
				# 1 : C will always be diagonal
				# >1: C is diagonal for first iterations, set to 0 afterwards
				if(flgDiagonalOnly < 1) flgDiagonalOnly = 0
				if(flgDiagonalOnly == 1){
					ccov1_sep = min(1, ccov1 * (N+1.5) / 3)
					ccovmu_sep = min(1-ccov1_sep, ccovmu * (N+1.5) / 3)
				} else if(N > 98 && flgdisplay && countiter == 0){
					cat('consider option DiagonalOnly for high-dimensional problems')
				}
				# ||ps|| is close to sqrt(mueff/N) for mueff large on linear fitness
				#damps = ... % damping for step size control, usually close to one 
				#    (1 + 2*max(0,sqrt((mueff-1)/(N+1))-1)) ... % limit sigma increase
				#    * max(0.3, ... % reduce damps, if max. iteration number is small
				#          1 - N/min(stopMaxIter,stopMaxFunEvals/lambda)) + cs; 
				damps = .eval.control('damps', cma)
				#qqq hacking of a different parameter setting, e.g. for ccov or damps,
				#  can be done here, but is not necessary anymore, see opts.CMA. 
				# ccov1 = 0.0*ccov1; disp(['CAVE: ccov1=' num2str(ccov1)]);
				# ccovmu = 0.0*ccovmu; disp(['CAVE: ccovmu=' num2str(ccovmu)]);
				# damps = inf*damps; disp(['CAVE: damps=' num2str(damps)]);
				# cc = 1; disp(['CAVE: cc=' num2str(cc)]);
			}
			########################
			# CONTINUE DEBUG HERE (798 matlab)
			# Display initial message
			if(countiter == 0 && flgdisplay){
				if(mu == 1){
					strw = '100'
				} else if(mu < 8){
					strw = c(sprintf("%.0f",100*weights[1]), "  ", sprintf("%.0f",100*weights[-1]))
				} else{
					strw = c(sprintf("%.2g",100*weights[1:3]), "...", tail(sprintf("%.2g",100*weights), 3))
				}
				cat("\nrun[",irun,"]  n=", N, ": (", mu, ",", lambda, ")-CMA-ES(w=", "[", strw,"], mu_eff=", round(mueff,2),")\n")
				if(flgDiagonalOnly == 1){
					cat("\n    C is diagonal")
				} else if(flgDiagonalOnly){
					cat(paste("\n    C is diagonal for ",  floor(flgDiagonalOnly), " iterations", sep = ""))
				}
			}
			countiter = countiter + 1
			# Generate and evaluate lambda offspring
			noiseEpsilon = noiseReevals = 0
			fitness$raw = rep(NA, lambda + noiseReevals)
			# parallel evaluation
			if(flgEvalParallel){
				arz = matrix(rnorm(N*lambda), N, lambda)
				if(flgDiagonalOnly!=1){
					arx =  repmat(as.matrix(xmean), 1, lambda) + sigma * (BD %*% arz)
					# Eq. (1)
				} else{
					# check
					arx = repmat(as.matrix(xmean), 1, lambda) + repmat(as.matrix(sigma * diagD), 1, lambda) * arz
				}
				#if(noiseHandling){
				#	if(noiseEpsilon == 0){
				#		arx = cbind(arx, arx[,1:noiseReevals])
				#	} else if(flgDiagonalOnly==1){
				#		arx = cbind(arx, arx[,1:noiseReevals] + repmat(noiseEpsilon * sigma * diagD, 1, noiseReevals) * matrix(rnorm(N*noiseReevals), N, noiseReevals))
				#	} else{
				#		arx = cbind(arx, arx[,1:noiseReevals] + noiseEpsilon * sigma * (BD %*% matrix(rnorm(N*noiseReevals), N, noiseReevals)))
				#	}
				#}
				# You may handle constraints here. You may either resample
				# arz(:,k) and/or multiply it with a factor between -1 and 1
				# (the latter will decrease the overall step size) and
				# recalculate arx accordingly. Do not change arx or arz in any
				# other way.
				if(!bnd$isactive){
					arxvalid = arx
				} else{
					arxvalid = xintobounds(arx, lbounds, ubounds)$x
				}
				# You may handle constraints here.  You may copy and alter
				# (columns of) arxvalid(:,k) only for the evaluation of the
				# fitness function. arx and arxvalid should not be changed.
				fitness$raw = fun(arxvalid, ...)
				#fitness$raw = fun(arxvalid, Data, forecast) 
				countevalNA = countevalNA + sum(which(is.na(fitness$raw)))
				counteval = counteval + length(na.omit(fitness$raw))
			}
			# non-parallel evaluation and remaining NaN-values
			# set also the reevaluated solution to NaN
			#if(noiseReevals>0) fitness$raw[lambda + which(is.na(fitness$raw(1:noiseReevals)))] = NA
			hn = which(is.na(fitness$raw))
			if(length(hn)>0){
				# AG: This is for the first run
				if(is.null(arxvalid)) arxvalid = matrix(NA, ncol = length(hn), nrow = N)
				if(is.null(arz)) arz = matrix(NA, ncol = length(hn), nrow = N)
				if(is.null(arx)) arx = matrix(NA, ncol = length(hn), nrow = N)
				# AG: This is the augmentation for subsequent re-runs
				if(dim(arx)[2]<length(hn)){
					diffl = length(hn) - dim(arx)[2]
					arx = cbind(arx, matrix(NA, ncol = diffl, nrow = N))
					arz = cbind(arz, matrix(NA, ncol = diffl, nrow = N))
					arxvalid = cbind(arxvalid, matrix(NA, ncol = diffl, nrow = N))
				}
				for(k in hn){
				# fitness.raw(k) = NaN; 
					tries = 0
					# Resample, until fitness is not NaN
					while(is.na(fitness$raw[k])){
						if(k <= lambda){
							# regular samples (not the re-evaluation-samples)
							arz[,k] = rnorm(N)
							# (re)sample
							if(flgDiagonalOnly){
								arx[,k] = xmean + sigma * diagD * arz[,k]
								# Eq. (1)
							} else{
								arx[,k] = xmean + sigma * (BD %*% arz[,k])
								# Eq. (1)
							}
						} else{
							# re-evaluation solution with index > lambda
							if(flgDiagonalOnly){
								arx[,k] = arx[,k-lambda] + (noiseEpsilon * sigma) * diagD * rnorm(N)
							} else{
								arx[,k] = arx[,k-lambda] + (noiseEpsilon * sigma) * (BD %*% rnorm(N))
							}
						}
						
						# You may handle constraints here. You may either resample
						# arz(:,k) and/or multiply it with a factor between -1 and 1
						# (the latter will decrease the overall step size) and
						# recalculate arx accordingly. Do not change arx or arz in any
						# other way.
						if(!bnd$isactive){
							arxvalid[,k] = arx[,k]
						} else{
							arxvalid[,k] = xintobounds(arx[,k], lbounds, ubounds)$x
						}
						# You may handle constraints here.  You may copy and alter
						# (columns of) arxvalid(:,k) only for the evaluation of the
						# fitness function. arx should not be changed.
						fitness$raw[k] = fun(arxvalid[,k], ...)
						#fitness$raw[k] = fun(arxvalid[,k], Data, forecast)
						tries = tries + 1
						if(is.na(fitness$raw[k])){
							countevalNA = countevalNA + 1
						}
						if(tries%%100 == 0){
							if(Warnings) warning(paste("\n", tries," NA objective function values at evaluation ", counteval, sep = ""))
						}
					}
					counteval = counteval + 1
				# retries due to NaN are not counted
				}
			}
			fitness$sel = fitness$raw
			#---- handle boundaries -----
			if(1 < 3 && bnd$isactive){
				# Get delta fitness values
				val = myprctile(fitness$raw, c(25, 75))
				# more precise would be exp(mean(log(diagC)))
				val = (val[2] - val[1]) / N / mean(diagC) / sigma^2
				#val = (myprctile(fitness.raw, 75) - myprctile(fitness.raw, 25)) ...
				#    / N / mean(diagC) / sigma^2;
				# Catch non-sensible values 
				if(!is.finite(val)){
					if(Warnings) warning("\nNon-finite fitness range")
					val = max(bnd$dfithist)
				} else if(val == 0){
					# happens if all points are out of bounds
					val = min(bnd$dfithist[which(bnd$dfithist>0)]) 
					# seems not to make sense, given all solutions are out of bounds
				} else if(bnd$validfitval == 0){
					# flag that first sensible val was found
					bnd$dfithist = NULL
					bnd$validfitval = 1
				}
				# Store delta fitness values
				if(length(bnd$dfithist) < (20+(3*N)/lambda)){
					bnd$dfithist = c(bnd$dfithist, val)
				} else{
					bnd$dfithist = c(bnd$dfithist[-1],val)
				}
				tmp = xintobounds(xmean, lbounds, ubounds)
				tx = tmp$x
				ti = tmp$idx
				# Set initial weights
				if(bnd$iniphase){
					if(length(ti)>0 && any(ti)!=0){
						bnd$weights[which(bnd$isbounded==0)] = 2.0002 * median(bnd$dfithist)
						if(bnd$flgscale == 0){
							#scale only initial weights then
							dd = diagC
							idx = which(bnd$isbounded==0)
							if(length(idx)>0){
								dd = dd[idx] / mean(dd)
								#  remove mean scaling
								bnd$weights[idx] = bnd$weights[idx]/ dd
							}
						}
						if(bnd$validfitval && countiter > 2){
							bnd$iniphase = 0
						}
					}
				}
				# Increase weights
				if(1 < 3 && length(ti)>0 && any(ti)!=0){
					# any coordinate of xmean out of bounds
					# judge distance of xmean to boundary
					tx = xmean - tx
					idx = (ti != 0 & abs(tx) > (3*max(1,sqrt(N)/mueff)) * sigma*sqrt(diagC))
					# only increase if xmean is moving away
					idx = idx & (sign(tx) == sign(xmean - xold))
					if(length(idx)>0){
						# increase
						# the factor became 1.2 instead of 1.1, because
						# changed from max to min in version 3.52
						bnd$weights[idx] = 1.2^(min(1, mueff/10/N)) * bnd$weights[idx]
					}
				}
				# Calculate scaling biased to unity, product is one
				if(bnd$flgscale != 0){ 
					bnd$scale = exp(0.9*(log(diagC)-mean(log(diagC))))
				}
				# Assigned penalized fitness
				bnd$arpenalty = as.numeric( (bnd$weights/bnd$scale) %*% (arxvalid - arx)^2)
				fitness$sel = fitness$raw + bnd$arpenalty
			}
			# handle boundaries
			# ----- end handle boundaries -----
			# Sort by fitness
			tmp = sort(fitness$raw, index.return = TRUE)
			fitness$raw = tmp$x
			fitness$idx = tmp$ix
			tmp = sort.int(fitness$sel, index.return = TRUE)
			fitness$sel = tmp$x
			fitness$idxsel = tmp$ix
			# minimization
			fitness$hist[-1] = fitness$hist[1:(length(fitness$hist)-1)]
			# record short history of
			fitness$hist[1] = fitness$raw[1]
			# best fitness values
			if(length(fitness$histbest) < (120+ceiling(30*N/lambda)) || (countiter%%5 == 0  && length(fitness$histbest) < 2e4)){
				# 20 percent of 1e5 gen.
				fitness$histbest = c(fitness$raw[1], fitness$histbest)
				# best fitness values
				fitness$histmedian = c(median(fitness$raw), fitness$histmedian)
				# median fitness values
			} else{
				fitness$histbest[-1] = fitness$histbest[1:(length(fitness$histbest)-1)]
				fitness$histmedian[-1] = fitness$histmedian[1:(length(fitness$histmedian)-1)]
				fitness$histbest[1] = fitness$raw[1]
				# best fitness values
				fitness$histmedian[1] = median(fitness$raw)
				# median fitness values
			}			
			fitness$histsel[-1] = fitness$histsel[1:(length(fitness$histsel)-1)]
			# record short history of best sel fitness values
			fitness$histsel[1] = fitness$sel[1]
			# Calculate new xmean, this is selection and recombination 
			xold = xmean
			# for speed up of Eq. (2) and (3)
			xmean = as.numeric(arx[,fitness$idxsel[1:mu]] %*% weights)
			zmean = as.numeric(arz[,fitness$idxsel[1:mu]] %*% weights)
			#==D^-1*B'*(xmean-xold)/sigma
			if(mu == 1){
				fmean = fitness$sel[1]
			} else{
				fmean = NA
				# [] does not work in the latter assignment
				# fmean = feval(fitfun, xintobounds(xmean, lbounds, ubounds), varargin{:});
				# counteval = counteval + 1;
			}	
			# Cumulation: update evolution paths
			ps = as.numeric( (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * (B%*%zmean) )
			# Eq. (4)
			hsig = as.integer(norm(as.matrix(ps), type = "F")/sqrt(1-(1-cs)^(2*countiter))/chiN < (1.4 + 2/(N+1)))
			#if(flg_future_setting){
			#	hsig = sum(ps^2) / (1-(1-cs)^(2*countiter)) / N < 2 + 4/(N+1)
			#	# just simplified
			#}
			#  hsig = norm(ps)/sqrt(1-(1-cs)^(2*countiter))/chiN < 1.4 + 2/(N+1);
			#  hsig = norm(ps)/sqrt(1-(1-cs)^(2*countiter))/chiN < 1.5 + 1/(N-0.5);
			#  hsig = norm(ps) < 1.5 * sqrt(N);
			#  hsig = 1;	
			pc = (1 - cc)*pc + hsig*(sqrt(cc*(2-cc)*mueff)/sigma) * (xmean-xold) 
			# Eq. (2)
			#if(hsig == 0){
			#	# disp([num2str(countiter) ' ' num2str(counteval) ' pc update stalled']);
			#}
			# Adapt covariance matrix
			neg$ccov = 0
			# TODO: move parameter setting upwards at some point
			if( (ccov1 + ccovmu) > 0){
				# Eq. (3)
				if(flgDiagonalOnly){
					#internal linear(?) complexity
					# regard old matrix
					# % plus rank one update
					#  % plus rank mu update
					diagC = (1-ccov1_sep-ccovmu_sep+(1-hsig)*ccov1_sep*cc*(2-cc)) * diagC + 
							ccov1_sep * pc^2 + ccovmu_sep * (diagC * (arz[,fitness$idxsel[1:mu]]^2 %*% weights))
					#   * (repmat(diagC,1,mu) .* arz(:,fitness.idxsel(1:mu)).^2 * weights);
					diagD = sqrt(diagC)
					# % replaces eig(C)
				} else{
					arpos = (arx[,fitness$idxsel[1:mu]] - repmat(as.matrix(xold),1,mu))/ sigma
					# "active" CMA update: negative update, in case controlling pos. definiteness
					if(flgActiveCMA > 0){
						# set parameters
						neg$mu = mu
						neg$mueff = mueff
						if(flgActiveCMA > 1){
							#flat weights with mu=lambda/2
							neg$mu = floor(lambda/2)
							neg$mueff = neg$mu
						}
						# neg.mu = ceil(min([N, lambda/4, mueff]));  neg.mueff = mu; % i.e. neg.mu <= N 
						# Parameter study: in 3-D lambda=50,100, 10-D lambda=200,400, 30-D lambda=1000,2000 a 
						# three times larger neg.ccov does not work. 
						#   increasing all ccov rates three times does work (probably because of the factor (1-ccovmu))
						#   in 30-D to looks fine
	
						neg$ccov = (1 - ccovmu) * 0.25 * neg$mueff / ((N+2)^1.5 + 2*neg$mueff)
						neg$minresidualvariance = 0.66
						# keep at least 0.66 in all directions, small popsize are most critical
						neg$alphaold = 0.5
						# where to make up for the variance loss, 0.5 means no idea what to do
						# 1 is slightly more robust and gives a better "guaranty" for pos. def., 
						# but does it make sense from the learning perspective for large ccovmu? 
						neg$ccovfinal = neg$ccov	
						# prepare vectors, compute negative updating matrix Cneg and checking matrix Ccheck
						arzneg = arz[,fitness$idxsel[seq(lambda, lambda - neg$mu + 1, by=-1)]]
						# i-th longest becomes i-th shortest
						# TODO: this is not in compliance with the paper Hansen&Ros2010, 
						#       where simply arnorms = arnorms(end:-1:1) ?
						tmp = sort.int(sqrt(apply(arzneg^2, 2, "sum")), index.return = TRUE)
						arnorms = tmp$x
						idxnorms = tmp$ix
						idxnorms = sort.int(idxnorms,  index.return = TRUE)$ix
						arnormfacs = arnorms[seq(length(arnorms), 1, by = -1)] / arnorms
						# arnormfacs = arnorms(randperm(neg.mu)) ./ arnorms;
						arnorms = arnorms[seq(length(arnorms), 1, by = -1)]
						# for the record
						if(flgActiveCMA < 2){
							arzneg = arzneg * repmat(t(arnormfacs[idxnorms]), N, 1)
							#  E x*x' is N
							# arzneg = sqrt(N) * arzneg ./ repmat(sqrt(sum(arzneg.^2, 1)), N, 1);  % E x*x' is N
						}
						if(flgActiveCMA < 1 && neg$mu == mu){
							# weighted sum
							if(flgActiveCMA%%1 == 1){
								#TODO: prevent this with a less tight but more efficient check (see below) 
								Ccheck = arzneg %*% diag(weights) %*% t(arzneg)
								# in order to check the largest EV
							}
							artmp = BD %*% arzneg
							Cneg = artmp %*% diag(weights) %*% t(artmp)
						} else{
							#simple sum
							if(flgActiveCMA%%1 == 1){
								Ccheck = (1/neg$mu) * arzneg%*%t(arzneg)
								# in order to check largest EV
							}
							artmp = BD %*% arzneg;
							Cneg = (1/neg$mu) * artmp%*%t(artmp)
						}
						# check pos.def. and set learning rate neg.ccov accordingly, 
						# this check makes the original choice of neg.ccov extremly failsafe 
						# still assuming C == BD*BD', which is only approxim. correct 
						if(flgActiveCMA%%1 == 1 && (1 - neg$ccov * arnorms[idxnorms]^2 %*% weights) < neg$minresidualvariance){
							# TODO: the simple and cheap way would be to set
							#    fac = 1 - ccovmu - ccov1 OR 1 - mueff/lambda and
							#    neg.ccov = fac*(1 - neg.minresidualvariance) / (arnorms(idxnorms).^2 * weights)
							# this is the more sophisticated way: 
							# maxeigenval = eigs(arzneg * arzneg', 1, 'lm', eigsopts);  # not faster
							maxeigenval = max(eigen(Ccheck)$values) 
							# norm is much slower, because (norm()==max(svd())
							#disp([countiter log10([neg.ccov, maxeigenval, arnorms(idxnorms).^2 * weights, max(arnorms)^2]), ...
							#          neg.ccov * arnorms(idxnorms).^2 * weights])
							# pause
							# remove less than ??34*(1-cmu)#?? of variance in any direction
							#     1-ccovmu is the variance left from the old C
							neg$ccovfinal = min(neg$ccov, (1-ccovmu)*(1-neg$minresidualvariance)/maxeigenval)
							# -ccov1 removed to avoid error message??
							if(neg$ccovfinal < neg$ccov){
								if(flgdisplay) cat(paste("\nactive CMA at iteration ",  countiter, " : max EV == ", sprintf("%0.3e", maxeigenval), " ", sprintf("%0.3e", neg$ccov), " ", sprintf("%0.3e",neg$ccovfinal), sep = ""))
							}
						}
						# xmean = xold;  # the distribution does not degenerate!? 
						# update C
						C = (1-ccov1-ccovmu+neg$alphaold*neg$ccovfinal+(1-hsig)*ccov1*cc*(2-cc)) * C + ccov1 *  pc%*%t(pc) + 
						(ccovmu + (1-neg$alphaold)*neg$ccovfinal) * arpos %*%  (repmat(as.matrix(weights),1,N) * t(arpos)) - neg$ccovfinal *  Cneg
						# minus rank mu update
					} else{  # no active (negative) update
						C = (1-ccov1-ccovmu+(1-hsig)*ccov1*cc*(2-cc)) * C + ccov1 *  pc%*%t(pc) + ccovmu * arpos %*% (repmat(as.matrix(weights),1,N) * t(arpos))
					# is now O(mu*N^2 + mu*N), was O(mu*N^2 + mu^2*N) when using diag(weights)
					#   for mu=30*N it is now 10 times faster, overall 3 times faster
					}
					diagC = diag(C)
				}
			}
			# the following is de-preciated and will be removed in future
			# better setting for cc makes this hack obsolete
			if(11 < 2 && !flgscience){
				# remove momentum in ps, if ps is large and fitness is getting worse.
				# this should rarely happen. 
				# this might very well be counterproductive in dynamic environments
				if(sum(ps^2)/N > (1.5 + 10*(2/N)^5) && fitness$histsel[1] > max(fitness$histsel[2:3])){
					ps = ps * sqrt(N*(1+max(0,log(sum(ps^2)/N))) / sum(ps^2))
					if(flgdisplay){
						cat(paste("\nMomentum in ps removed at [niter neval]=[", countiter, " ", counteval,"]",sep=""))
					}
				}
			}
			# Adapt sigma
			if(flg_future_setting){
				# according to a suggestion from Dirk Arnold (2000)
				# exp(1) is still not reasonably small enough
				sigma = sigma * exp(min(1, (sum(ps^2)/N - 1)/2 * cs/damps))
				# Eq. (5)
			} else{
				# exp(1) is still not reasonably small enough
				sigma = sigma * exp(min(1, (sqrt(sum(ps^2))/chiN - 1) * cs/damps))
				# Eq. (5)
			}
			# disp([countiter norm(ps)/chiN])
			if(11 < 3){
				# testing with optimal step-size
				if(countiter == 1){
					if(flgdisplay) cat('*********** sigma set to const * ||x|| ******************')
				}
				sigma = 0.04 * mueff * sqrt(sum(xmean^2))/N
				# 20D,lam=1000:25e3
				sigma = 0.3 * mueff * sqrt(sum(xmean^2))/N
				# 20D,lam=(40,1000):17e3
				#      75e3 with def (1.5)
				#      35e3 with damps=0.25
			}
			if(11 < 3){
				if(countiter == 1){
					if(flgdisplay) cat('\n*********** xmean set to const ******************')
				}
				xmean = ones(N,1)
			}
			# Update B and D from C
			
			if(!flgDiagonalOnly && (ccov1+ccovmu+neg$ccov) > 0 && (countiter%%(1/(ccov1+ccovmu+neg$ccov)/N/10)) < 1){
				C = triu(C) + t(triu(C,1))
				# enforce symmetry to prevent complex numbers
				xtmp = eigen(C)
				B = xtmp$vector
				tmp = xtmp$values
				# eigen decomposition, B==normalized eigenvectors
				# effort: approx. 15*N matrix-vector multiplications
				diagD = tmp
				if(any(!is.finite(diagD))){
					stop(paste("\nfunction eigen returned non-finited eigenvalues, cond(C)=", cond(C), sep = ""))
				}
				if(any(!is.finite(as.vector(B)))){
					stop(paste("\nfunction eigen returned non-finited eigenvectors, cond(C)=", cond(C), sep = ""))
				}
				# limit condition of C to 1e14 + 1
				if(min(diagD) <= 0){
					if(stopOnWarnings){
						stopflag = c(stopflag, "warnconditioncov")
					} else{
						if(Warnings) warning(paste("\nIteration ", countiter," :Eigenvalue (smaller) than zero", sep = ""))
						diagD[diagD<0] = 0
						tmp = max(diagD)/1e14
						C = C + tmp*eye(N,N)
						diagD = diagD + tmp*ones(N,1)
					}
				}
				if(max(diagD) > 1e14*min(diagD)){
					if(stopOnWarnings){
						stopflag = c(stopflag, "warnconditioncov")
					} else{
						if(Warnings) warning(paste("\nIteration ", countiter," :condition of C at upper limit", sep = ""))
						tmp = max(diagD)/1e14 - min(diagD)
						C = C + tmp * eye(N,N)
						diagD = diagD + tmp*ones(N,1)
					}
				}
				diagC = diag(C)
				diagD = sqrt(diagD)
				# D contains standard deviations now
				# diagD = diagD / prod(diagD)^(1/N);  C = C / prod(diagD)^(2/N);
				BD = B*repmat(t(diagD),N,1)
				# O(n^2)
			}		
			# Align/rescale order of magnitude of scales of sigma and C for nicer output
			# not a very usual case
			if(1 < 2 && sigma > 1e10*max(diagD)){
				fac = sigma / max(diagD)
				sigma = sigma/fac
				pc = fac * pc
				diagD = fac * diagD
				if(!flgDiagonalOnly){
					C = fac^2 * C
					# disp(fac);
					BD = B*repmat(t(diagD),N,1)
					# # O(n^2), but repmat might be inefficient todo?
				}
				diagC = fac^2 * diagC
			}
			if(flgDiagonalOnly > 1 && countiter > flgDiagonalOnly){
				# full covariance matrix from now on 
				flgDiagonalOnly = 0
				B = eye(N,N)
				BD = diag(diagD)
				C = diag(diagC)
				# is better, because correlations are spurious anyway
			}							
									
			# ----- numerical error management -----
			# Adjust maximal coordinate axis deviations
			if(any(sigma*sqrt(diagC) > maxdx)){
				sigma = min(maxdx /sqrt(diagC))
			}
			# Adjust minimal coordinate axis deviations
			if(any(sigma*sqrt(diagC) < mindx)){
				sigma = max(mindx/sqrt(diagC)) * exp(0.05+cs/damps)
			}
			# Adjust too low coordinate axis deviations
			if(any(xmean == (xmean + 0.2*sigma*sqrt(diagC)))){
				if(stopOnWarnings){
					stopflag = c(stopflag, "'warnnoeffectcoord")
				} else{
					if(Warnings) warning(paste("\nIteration ", countiter, ": coordinate axis std deviation too low", sep=""))
					if(flgDiagonalOnly){
						diagC = diagC + (ccov1_sep+ccovmu_sep) * (diagC * (xmean == xmean + 0.2*sigma*sqrt(diagC)))
					} else{
						C = C + (ccov1+ccovmu) * diag(diagC * (xmean == xmean + 0.2*sigma*sqrt(diagC)))
					}
				sigma = sigma * exp(0.05+cs/damps)
				}
			}
			# Adjust step size in case of (numerical) precision problem 
			if(flgDiagonalOnly){
				tmp = 0.1*sigma*diagD
			} else{
				tmp = 0.1*sigma*BD[,1+floor(countiter%%N)]
			}
			if(all(xmean == (xmean + tmp))){
				i = 1+floor(countiter%%N)
				if(stopOnWarnings){
					stopflag = stopflag(stopflag, "warnnoeffectaxis")
				} else{
					warning(paste("\nIteration ", countiter, ": main axis standard deviation ", sigma*diagD[i], " has no effect", sep = ""))
					sigma = sigma * exp(0.2+cs/damps)
				}
			}
			# Adjust step size in case of equal function values (flat fitness)
			# isequalfuncvalues = 0; 
			if(fitness$sel[1] == fitness$sel[1+ceiling(0.1+lambda/4)]){
				# isequalfuncvalues = 1; 
				if(stopOnEqualFunctionValues==1){
					nx = length(as.numeric(arrEqualFunvals))
					arrEqualFunvals = c(countiter, arrEqualFunvals[1:(nx-1)])
					# stop if this happens in more than 33#
					if(arrEqualFunvals[nx] > (countiter - 3 * nx)){
						stopflag = c(stopflag, "equalfunvals")
					}
				} else{
					if(flgWarnOnEqualFunctionValues){
						if(Warnings) warning(paste("\nIteration ", countiter, ": equal function values f=",  fitness$sel[1], " at maximal main axis sigma ", sigma*max(diagD), sep = ""))
					}
					sigma = sigma * exp(0.2+cs/damps)
				}
			}
			# Adjust step size in case of equal function values
			
			if(countiter > 2 && myrange(c(fitness$hist, fitness$sel[1])) == 0){
				if(stopOnWarnings){
					stopflag = c(stopflag, "warnequalfunvalhist")
				} else{
					if(Warnings) warning(paste("\nIteration ", countiter, ": equal function values in history at maximal main axis sigma ", sigma*max(diagD), sep = ""))
				}
				sigma = sigma * exp(0.2+cs/damps)
			}
			# ----- end numerical error management ----
			# Keep overall best solution
			out$evals = counteval
			out$solutions$evals = counteval
			out$solutions$mean$x = xmean
			out$solutions$mean$f = fmean
			out$solutions$mean$evals = counteval
			out$solutions$recentbest$x = arxvalid[,fitness$idx[1]]
			out$solutions$recentbest$f = fitness$raw[1]
			out$solutions$recentbest$evals = counteval + fitness$idx[1] - lambda
			out$solutions$recentworst$x = arxvalid[, fitness$idx[length(fitness$idx)]]
			out$solutions$recentworst$f = fitness$raw[length(fitness$raw)]
			out$solutions$recentworst$evals = counteval + fitness$idx[length(fitness$idx)] - lambda
			if(fitness$hist[1] < out$solutions$bestever$f){
				out$solutions$bestever$x = arxvalid[,fitness$idx[1]]
				out$solutions$bestever$f = fitness$hist[1]
				out$solutions$bestever$evals = counteval + fitness$idx[1] - lambda
				bestever = out$solutions$bestever
			}
			# Set stop flag
			if(fitness$raw[1] <= stopFitness){
				stopflag = c(stopflag, "fitness")
			}
			if(counteval >= stopMaxFunEvals){
				stopflag = c(stopflag, "maxfunevals")
			}	
			if(countiter >= stopMaxIter){
				stopflag = c(stopflag, "maxiter")
			}	
			if(all(sigma*(max(abs(pc), sqrt(diagC))) < stopTolX)){
				stopflag = c(stopflag, "tolx")
			}
			if(any(sigma*sqrt(diagC) > stopTolUpX)){
				stopflag = c(stopflag, "tolupx")
			}
			if(sigma*max(diagD) == 0){
				# should never happen
				stopflag = c(stopflag, "bug")
			}
			if(countiter > 2 && myrange(c(fitness$sel, fitness$hist)) <= stopTolFun){
				stopflag = c(stopflag, "tolfun")
			}
			if(countiter >= length(fitness$hist) && myrange(fitness$hist) <= stopTolHistFun){
				stopflag = c(stopflag, "tolhistfun")
			}
			l = floor(length(fitness$histbest)/3)
			nl = length(fitness$histmedian)
			if(1 < 2 && stopOnStagnation && countiter > (N * (5+100/lambda)) && length(fitness$histbest) > 100 && 
					median(fitness$histmedian[1:l]) >= median(fitness$histmedian[(nl-l):nl]) && 
					median(fitness$histbest[1:l]) >= median(fitness$histbest[(nl-l):nl])){
				stopflag = c(stopflag, "stagnation")		
			}
			if(counteval >= stopFunEvals || countiter >= stopIter){
				stopflag = c(stopflag, "stoptoresume")		
				if(length(stopflag) == 1 && flgsaving == 0){
					stop('To resume later the saving option needs to be set')
				}
			}
		
			out$stopflag = stopflag
			
			# ----- output generation -----
			if(verbosemodulo > 0 && is.finite(verbosemodulo)){
				if(countiter == 1 || countiter%%(10*verbosemodulo) < 1){ 
					#cat("\nIter, #f:   f-value\t  (median,worst)\t |Axis Ratio|\t idx:Min SD idx:Max SD")
					cat("\n Iter        #f:           f-value")
				}
				if(countiter%%verbosemodulo < 1 || (verbosemodulo > 0 && is.finite(verbosemodulo) && 
					(countiter < 3 || !is.null(stopflag)))){
					tmp = min.int(sigma*sqrt(diagC), index.return = TRUE)
					minstd  = tmp$x
					minstdidx = tmp$ix
					tmp = max.int(sigma*sqrt(diagC), index.return = TRUE)
					maxstd  = tmp$x
					maxstdidx = tmp$ix			
					# format display nicely
					#cat(paste("\n",rep(" ", 4-floor(log10(countiter))), countiter, rep(" ", 1,5-floor(log10(counteval))), 
					#		",",counteval," : ", sprintf('%.5e', fitness$hist[1]), " + (", sprintf('%.0e', median(fitness$raw)-fitness$hist[1],4)," , ", 
					#				sprintf('%.0e', max(fitness$raw)-fitness$hist[1],4), ")  | ",  sprintf('%4.2e', max(diagD)/min(diagD),4), " |   ", 
					#				rep(" ",1-floor(log10(minstdidx))), minstdidx, ":", sprintf('%.1e,',minstd,4), " ",
					#				rep(" ",1-floor(log10(maxstdidx))), maxstdidx, ":", sprintf('%.1e,', maxstd,4), sep = ""))
					cat(c("\n",countiter, rep("", max(1, 10-floor(log10(countiter)))), counteval, 
									rep("", max(1, 10-floor(log10(counteval)))), sprintf('%.5e', fitness$hist[1])))
				}
				# measure time for recording data
				if(countiter < 3){
					timer$c = 0.05
					timer$nonoutput = 0
					timer$recording = 0
					timer$saving  = 0.15
					# first saving after 3 seconds of 100 iterations
					timer$plotting = 0
				} else if(countiter > 300){
					# set backward horizon, must be long enough to cover infrequent plotting etc
					# time.c = min(1, time.nonoutput/3 + 1e-9); 
					timer$c = max(1e-5, 0.1/sqrt(countiter))
					# mean over all or 1e-5
				}
				# get average time per iteration
				timer$t1 = Sys.time()
				timer$act = max(0, timer$t1 - timer$t0)
				timer$nonoutput = (1-timer$c) * timer$nonoutput + timer$c * timer$act
				timer$recording = (1-timer$c) * timer$recording
				# per iteration
				timer$saving   = (1-timer$c) * timer$saving
				timer$plotting = (1-timer$c) * timer$plotting	
			
				# get average time for recording data
				timer$t2 = Sys.time()
				timer$recording = timer$recording + timer$c * max(0, timer$t2 - timer$t1) 
				# plot
				timer$t3 = Sys.time()
				if(!is.null(stopflag) || timer$saving < (0.05 * timer$nonoutput) || countiter == 100){
					xmin = arxvalid[, fitness$idx[1]]
					fmin = fitness$raw[1]
					timer$saving = timer$saving + timer$c * max(0,Sys.time() - timer$t3) 
					}
				}
				timer$t0 = Sys.time()
			}
		}
		fmin = fitness$raw[1]
		xmin = arxvalid[, fitness$idx[1]]
		# Return best point of last generation.
		dix = which(stopflag == "stoptoresume")
		if(length(dix)>0 && length(stopflag) > length(dix)){
			# final stopping
			out$solutions$mean$f = fun(xintobounds(xmean, lbounds, ubounds), ...)
			counteval = counteval + 1
			out$solutions$mean$evals = counteval
			if(out$solutions$mean$f < fitness$raw[1]){
				fmin = out$solutions$mean$f
				xmin = xintobounds(xmean, lbounds, ubounds)
				#% Return xmean as best point
			}
			if(out$solutions$mean$f < out$solutions$bestever$f){
				out$solutions$bestever = out$solutions$mean
				# Return xmean as bestever point
				out$solutions$bestever$x = xintobounds(xmean, lbounds, ubounds)
				bestever = out$solutions$bestever
			}
		}
		cat("\n")
		return(list(par = xmin, objective = fmin, out = out, counteval = counteval, stopflag = stopflag, bestever = bestever))
}

################################################################################
# Test Functions
frosenbrock = function(x){
	f = 100*sum((x[1:(length(x)-1)]^2 - x[-1])^2) + sum((x[1:(length(x)-1)]-1)^2)
	return(f)
}

fsphere = function(x){ return( sum(x^2) ) }

fssphere = function(x){ return( sqrt(sum(x^2)) ) }

fschwefel = function(x){
	f = 0
	for(i in 1:length(x)){
		f = f+sum(x[1:i])^2
	}
	return(f)
}

fcigar = function(x){ return( x[1]^2 + 1e6*sum(x[-1]^2) ) }

fcigtab = function(x){ return( x[1]^2 + 1e8*x[length(x)]^2 + 1e4*sum(x[2:(length(x)-1)]^2) ) }

ftablet = function(x){ return( 1e6*x[1]^2 + sum(x[-1]^2) ) }

felli = function(x){
	N = length(x)
	f = sum(1e6^((0:(N-1))/(N-1)) * x^2)
	return(f)
}

felli100 = function(x){
	N = length(x)
	f = sum(1e4^((0:(N-1))/(N-1)) * x^2)
	return(f)
}

fplane = function(x){ return(x[1]) }

ftwoaxes = function(x){
	N = length(x)
	if(N==1) stop("\nlength of x must be greater than 1")
	f = sum(x[1:floor(N/2)]^2) + 1e6*sum(x[floor(1+N/2):N]^2)
	return(f)
}

fparabR = function(x){
	f = -x[1] + 100*sum(x[-1]^2)
	return(f)
}

fsharpR = function(x){
	f = -x[1] + 100*norm(matrix(x[-1]), type = "f")
	return(f)
}

fdiffpow = function(x){
	N = length(x)
	f = sum(abs(x)^(2+10*(0:(N-1))/(N-1)))
	return(f)
}
					
frastrigin10 = function(x){
	N = length(x)
	scale = 10^((0:(N-1))/(N-1))
	f = 10*N + sum((scale*x)^2 - 10*cos(2*pi*(scale*x)))
	return(f)
}