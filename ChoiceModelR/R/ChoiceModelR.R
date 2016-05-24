choicemodelr <-
function(data, xcoding, demos, prior, mcmc, constraints, options) {
	callStop = function(message) { stop(message, call. = FALSE) }

	if (missing(data)) { callStop("data argument required") }
	if (ncol(data) < 5) { callStop("there are not enough columns in data to specify ID, Set, Alt, X variable(s), and y") }
	natts = ncol(data) - 4
	if (missing(xcoding)) { callStop("xcoding argument required") }
	if (length(xcoding) != natts) { callStop("length of xcoding is not equal to number of attributes") }
	if (any(xcoding != 0 & xcoding != 1)) { callStop("xcoding must contain only values 0 and 1") }
	npw = sum(xcoding == 0)
	ID = unique(data[,1])
	nunits = length(ID)
	drawdelta = ifelse(missing(demos), FALSE, TRUE)
	if (drawdelta) {
		if (nrow(demos) != nunits) { callStop("number of rows in demos does not equal number of units") }
		nz = ncol(demos)
		demos = t(t(demos) - colMeans(demos))
	}
	if (missing(options)) {
		none = FALSE
		save = FALSE
		restart = FALSE
		keep = 10
		wgt = 5
	}
	else {
		none = ifelse(is.null(options$none), FALSE, options$none)
		save = ifelse(is.null(options$save), FALSE, options$save)
		restart = ifelse(is.null(options$restart), FALSE, options$restart)
		keep = ifelse(is.null(options$keep), 10, options$keep)
		wgt = ifelse(is.null(options$wgt), 5, options$wgt)
		if (wgt < 1 | wgt > 10) { callStop("wgt must be in the range 1 to 10") }
	}

	maxsets = max(data[,2])
	maxalts = max(data[,3]) + ifelse(none, 1, 0)

	nsets = rep(0, nunits)
	nalts = matrix(0, nrow = nunits, ncol = maxsets)
	for (i in 1:nunits) {
		temp = data[data[,1] == ID[i],]
		nsets[i] = max(temp[,2])
		for (j in 1:nsets[i]) {
			nalts[i, j] = max(temp[temp[,2] == j, 3]) + ifelse(none, 1, 0)
		}
	}

	setsmap = (nalts > 0)
	altsmap = matrix(FALSE, nrow = sum(nsets), ncol = maxalts)
	index = 1
	nalts = t(nalts) # TO GET AROUND R COLUMN FILL
	for (i in 1:length(nalts)) {
		if (nalts[i] != 0) {
			altsmap[index, 1:nalts[i]] = TRUE
			index = index + 1
		}
	}

	avgalts = mean(nalts[nalts > 0])
	avgsets = mean(nsets)
	nalts = colSums(nalts)
	
	info = list(nunits = nunits,
			maxsets = maxsets,
			maxalts = maxalts,
			nsets = nsets,
			nalts = nalts,
			setsmap = setsmap,
			altsmap = altsmap)	

	if (all(data[!(data[,3] == 1), ncol(data)] == 0)) { share = FALSE }
	else { share = TRUE }

	if (share) {
		Xt = data[, c(-1, -2, -3, -ncol(data)), drop = FALSE]
		y = matrix(0, nrow = maxalts, ncol = sum(nsets))
		y[t(altsmap)] = data[, ncol(data)]
		y = t(y)
		ytab = colSums(y > 0); ytab = rbind(ytab, round(ytab / sum(ytab) * 100, 2))
		y = y / rowSums(y) * wgt
	}
	else {
		if (none) {
			Xt = NULL
			for (i in 1:nunits) {
				for (j in 1:nsets[i]) {
					ind = (data[,1] == ID[i]) & (data[,2] == j)
					Xt = rbind(Xt, cbind(data[ind, c(-1, -2, -3, -ncol(data)), drop = FALSE], 0), c(rep(0, natts), 1))
				}
			}
		}
		else { Xt = data[, c(-1, -2, -3, -ncol(data)), drop = FALSE] }
		y = data[data[,3] == 1, ncol(data)]
		if (any(y > maxalts) | any(y < 1)) { callStop(paste("invalid values of y present in data - values must be 1 to ", 
									    maxalts, sep = "")) }
		ytab = table(y); ytab = rbind(ytab, round(ytab / sum(ytab) * 100, 2))
	}
	
	nlevels = c(rep(0, natts))
	for (i in 1:natts) {
		if (xcoding[i] == 0) {
			nlevels[i] = length(unique(Xt[,i]))
			if (any(unique(Xt[,i]) == 0)) { nlevels[i] = nlevels[i] - 1 }
		}
		else { nlevels[i] = 1 }
	}
	npar = sum(nlevels) - npw + ifelse(none, 1, 0)
	X = matrix(0, nrow = nrow(Xt), ncol = npar)
	Xind = 1
	if (none) { nn = Xt[, ncol(Xt)] == 0 }
	else { nn = rep(TRUE, nrow(Xt)) }
	for (i in 1:natts) {
		if (xcoding[i] == 1) {
			X[nn, Xind] = Xt[nn, i] - mean(unique(Xt[nn, i]))
			Xind = Xind + 1
		}
		else {
			for (j in 1:(nlevels[i] - 1)) {
				X[nn, Xind] = (Xt[nn, i] == j) * 1 - (Xt[nn, i] == nlevels[i]) * 1
				Xind = Xind + 1
			}
		}
	}
	if (none) { X[, npar] = Xt[, ncol(Xt)] }

	effectsmap = matrix(FALSE, nrow = natts, ncol = npar)
	index = 1
	for (i in 1:natts) {
		count = 1
		repeat {
			effectsmap[i, index] = TRUE
			index = index + 1
			count = count + 1
			if (count > max(1, nlevels[i] - 1)) { break }
		}
	}
	if (none) { effectsmap = rbind(effectsmap, c(rep(FALSE, npar - 1), TRUE)) }

	if (missing(prior)) { prior = list() }
	if (is.null(prior$mubar)) { mubar = matrix(rep(0, npar), nrow = 1) }
	else {
		mubar = matrix(prior$mubar, nrow = 1)
		if (ncol(mubar) != npar) { callStop("ncol(mubar) does not equal npar") }
	}
	if (is.null(prior$Amu)) { Amu = matrix(0.01, ncol = 1) }
	else {
		Amu = matrix(prior$Amu, ncol = 1)
		if (any(dim(Amu) != c(1, 1))) { callStop("Amu must be a 1 x 1 matrix") }
	}
	if (is.null(prior$df)) { df = 5 }
	else {
		df = prior$df
		if (df < 2) { callStop("invalid value of df - must be >= 2") }
	}
	nu = npar + df
	if (is.null(prior$v)) {
		v = 2
	}
	else {
		v = prior$v
		if (v <= 0) { callStop("invalid v - must be >= 0") }
	}
	if (is.null(prior$Ad) & drawdelta) { Ad = 0.01 * diag(npar * nz) }
	else if (drawdelta) {
		Ad = prior$Ad
		if (any(dim(Ad) != c(npar * nz, npar * nz))) { callStop("Ad must be of dimensions (npar * nz) x (npar * nz)") }
	}
	if (is.null(prior$deltabar) & drawdelta) { deltabar = rep(0, nz * npar) }
	else if (drawdelta) {
		deltabar = prior$deltabar
		if (length(deltabar) != nz * npar) { callStop("deltabar must be of length npar * nz") }
	}

	V = matrix(0, nrow = npar, ncol = npar)
	v.ind = 1
	for (i in 1:natts) {
		if (nlevels[i] == 1) {
			V[v.ind, v.ind] = 1
			v.ind = v.ind + 1
		}
		else {
			pcov = -1 / nlevels[i]
			pvar = (nlevels[i] - 1) / nlevels[i]
			temp = pvar * diag(nlevels[i] - 1)
			temp[upper.tri(temp) | lower.tri(temp)] = pcov

			V[v.ind:(v.ind + nlevels[i] - 2), v.ind:(v.ind + nlevels[i] - 2)] = temp
			v.ind = v.ind + nlevels[i] - 1
		}
	}
	if (none) { V[v.ind, v.ind] = 1 }
	V = V * nu * v

	if (missing(mcmc)) { callStop("mcmc argument required") }
	if (is.null(mcmc$s)) {
		s = 0.1
		adjust.s = TRUE
	}
	else {
		s = mcmc$s
		adjust.s = FALSE
	}
	if (is.null(mcmc$R)) { callStop("element R of mcmc required") }
	R = mcmc$R
	if (is.null(mcmc$use)) { use = min(0.5 * R, 10000) }
	else {
		use = mcmc$use
		if (use > R) { callStop("use must be <= R") }
	}

	constrain = ifelse(missing(constraints), FALSE, TRUE)
	if (constrain) {
		if (length(constraints) != natts) { callStop("length(constraints) must equal number of attributes") }
		for (i in 1:natts) {
			if (any(dim(constraints[[i]]) != c(nlevels[i], nlevels[i]))) {
				callStop("dim(constraints[[i]]) must equal nlevels[i] x nlevels[i]")
			}
			if (any(constraints[[i]] != -1 & constraints[[i]] != 0 & constraints[[i]] != 1)) {
				callStop("constraints matrices must only contain -1, 0, and 1 elements")
			}
		}
	}

	fR = 0
	if (restart) {
		if (!any(list.files() == "restart.txt")) {
			callStop("file restart.txt does not exist in the working directory")
		}
		f.in = scan("restart.txt", nlines = 1, skip = 0, quiet = TRUE)
		fR = f.in[1]
		fUnits = f.in[2]
		fPar = f.in[3]
		s = f.in[4]
		fCons = f.in[5]
		fDem = f.in[6]
		items = 1
		ltr = fUnits
		if (fCons == 1) { 
			items = items + 1 
			ltr = c(ltr, fUnits)	
		} 
		if (fDem == 1) { 
			items = items + 1 
			ltr = c(ltr, 1)
		}
		if (fUnits != nunits) { callStop("number of units in restart file does not match data file") }	
		if (fPar != npar) { callStop("number of parameters in restart file does not match data file") }
		if (fCons != as.numeric(constrain)) { callStop("restart file and function arguments do not agree on constraints being present") }
		if (fDem != as.numeric(drawdelta)) { callStop("restart file and function arguments do not agree on demographics being present") }
	
		fIndex = 1
		fCount = 1
		repeat {
			if (fCount > items) { break }
			f.in = scan("restart.txt", nlines = ltr[fCount], skip = fIndex, quiet = TRUE)
			switch(fCount,
				lb <- matrix(f.in, ncol = fPar, byrow = TRUE),
				lbc <- matrix(f.in, ncol = fPar, byrow = TRUE),
				ld <- c(f.in)
			)
			fIndex = fIndex + ltr[fCount]
			fCount = fCount + 1
		}
	}

	if (save) {
		if (drawdelta) { deltadraw = matrix(0, nrow = floor(use / keep), ncol = nz * npar) }
		betadraw = array(0, dim = c(nunits, npar, floor(use / keep)))
		compdraw = NULL
		loglike = rep(0, floor(use / keep))
	}
	if (constrain) {
		if (save) { betadraw.c = array(0, dim = c(nunits, npar, floor(use / keep))) }
		if (restart) { oldbetas.c = lbc }
		else { oldbetas.c = matrix(0, nrow = nunits, ncol = npar) }
	}
	betaout = matrix(0, nrow = nunits, ncol = npar)
	if (drawdelta) { 
		if (restart) { olddelta = ld } 
		else { olddelta = rep(0, nz * npar) }
	}
	if (restart) { oldbetas = lb }
	else { oldbetas = matrix(0, nrow = nunits, ncol = npar) }
	oldll = rep(0, nunits)
	oldcomp = NULL
	
	muplot = matrix(0, nrow = R, ncol = npar + npw)
	xplot = (1 + fR):(R + fR)

	#
	# DEFINE NEEDED FUNCTIONS
	#
	getLLMnl = function (beta, y, X, info) {
		nunits = info$nunits
		nalts = info$maxalts
		nsets = info$maxsets

		map.sets = t(info$setsmap)
		map.alts = t(info$altsmap)
		
		tmpXbeta = X * beta[rep(1:nunits, info$nalts),]
		tmpXbeta = rowSums(tmpXbeta)
		Xbeta = matrix(0, nrow = nalts, ncol = sum(info$nsets))
		Xbeta[map.alts] = tmpXbeta

		# SHARE DATA
		if (is.matrix(y)) {
			Xbeta[map.alts] = exp(Xbeta[map.alts])
			denom = colSums(Xbeta)
			probs = (t(Xbeta) / denom) ^ y
			tmpProbs = apply(probs, 1, prod)
			probs = matrix(1, nrow = nsets, ncol = nunits)
			probs[map.sets] = tmpProbs

			ll = log(apply(probs, 2, prod))
		}

		# DISCRETE DATA
		else {
			ind = cbind(y, 1:sum(info$nsets))
			xby = matrix(0, nrow = nsets, ncol = nunits)
			xby[map.sets] = Xbeta[ind]
			Xbeta[map.alts] = exp(Xbeta[map.alts])
			denom = matrix(0, nrow = nsets, ncol = nunits)
			denom[map.sets] = log(colSums(Xbeta))
			
			ll = colSums(xby - denom)
		}

		return(ll)
	}

	getLndMvn = function (x, mu, rooti) {
		npar = ncol(x)
		#
		# WITH COVARIATES
		#	
		if (is.matrix(mu)) { z = (x - mu) %*% rooti }
		#
		# NO COVARIATES
		#
		else { z = crossprod(t(x) - mu, rooti) }

		logs = -(npar / 2) * log(2 * pi) - 0.5 * rowSums(z * z) + sum(log(diag(rooti)))
		
		return(logs)
	}

	drawDelta = function (x, y, comps, deltabar, Ad) {
		yy = t(t(y) - comps$mu)
		sig = tcrossprod(comps$rooti)

		xtx = crossprod(x) %x% sig
	        xty = matrix(sig %*% crossprod(yy, x), ncol = 1)

		cov = chol2inv(chol(xtx + Ad))

		return(cov %*% (xty + Ad %*% deltabar) + t(chol(cov)) %*% rnorm(length(deltabar)))
	}

	mnlRwMetropOnce = function (y, X, oldbeta, oldll, s, inc.root, betabar, rootpi, info, constraints, oldbeta.c) {
		stay = 0
		nunits = info$nunits
		npar = ncol(oldbeta)

		increment = s * crossprod(matrix(rnorm(nunits * npar), ncol = nunits), inc.root)
		newbeta = oldbeta + increment

		if (!missing(constraints)) {
			newbeta.c = constrainBetas(newbeta, constraints)
			newll = getLLMnl(newbeta.c, y, X, info)
		}
		else {
			newll = getLLMnl(newbeta, y, X, info)
		}

		newlpost = newll + getLndMvn(newbeta, betabar, rootpi)
		ldiff = newlpost - oldll - getLndMvn(oldbeta, betabar, rootpi)

		alpha = exp(ldiff)
		alpha[alpha > 1] = 1
		unif = runif(nunits)
		unif[alpha == 1] = 0
	
		good = (unif <= alpha)
		betadraw = oldbeta
		betadraw[good,] = newbeta[good,]
		if (!missing(constraints)) {
			betadraw.c = oldbeta.c
			betadraw.c[good,] = newbeta.c[good,]
		}
		oldll[good] = newll[good]
		
		stay = sum(!good)

		if (!missing(constraints)) {
			return(list(betadraw = betadraw, betadraw.c = betadraw.c, stay = stay, oldll = oldll))
		}
		else {
			return(list(betadraw = betadraw, stay = stay,  oldll = oldll))
		}
	}

	constrainBetas = function (betas, constraints) {
		newlength = 0
		oldlength = 0
		for (h in 1:length(constraints)) {
			if (all(constraints[[h]] == 0)) {
				if (ncol(constraints[[h]]) == 1) {
					newlength = 1
				}
				else {
					newlength = ncol(constraints[[h]]) - 1
				}
				oldlength = oldlength + newlength
				next
			}
			else if (ncol(constraints[[h]]) == 1) {
				newlength = 1
				betat = betas[,oldlength + 1]
				if (constraints[[h]] == 1) {
					betat[betat < 0] = 0
				}
				if (constraints[[h]] == -1) {
					betat[betat > 0] = 0
				}
				betas[,oldlength + 1] = betat
			}
			else {
				newlength = ncol(constraints[[h]]) - 1
				betat = betas[,(oldlength + 1):(oldlength + newlength), drop = FALSE]
				betat = cbind(betat, -1 * rowSums(betat))
				repeat {
					bad = FALSE
					for (i in 1:(nrow(constraints[[h]]) - 1)) {
						for (j in (i + 1):ncol(constraints[[h]])) {
							if (constraints[[h]][i, j] == 1) {
								change = (betat[,i] < betat[,j])
								if (sum(change) > 0) {
									betat[change, j] = betat[change, i]
									betat[change,] = betat[change,] - rowMeans(betat[change,,drop = FALSE])
									bad = TRUE
								}
							}
							else if (constraints[[h]][i, j] == -1) {
								change = (betat[,i] > betat[,j])
								if (sum(change) > 0) {
									betat[change, i] = betat[change, j]
									betat[change,] = betat[change,] - rowMeans(betat[change,,drop = FALSE])
									bad = TRUE
								}
							}
						}
					}
					if (!bad) { break }
				}
				betas[,(oldlength + 1):(oldlength + newlength)] = betat[,1:newlength]
			}
			oldlength = oldlength + newlength
		}
		return(betas)
	}

	rGibbs = function (y, betabar, A, nu, V) {
		temp = rmultireg(y, matrix(rep(1, nrow(y)), ncol = 1), betabar, A, nu, V)
		comps = list(mu = as.vector(temp$B), rooti = backsolve(chol(temp$Sigma), diag(ncol(temp$Sigma))))

		return(comps)
	}

	rmultireg = function (Y, X, Bbar, A, nu, V) {
   		n = nrow(Y)
    		m = ncol(Y)
	    	k = ncol(X)
    		RA = chol(A)
    		W = rbind(X, RA)
    		Z = rbind(Y, RA %*% Bbar)
    		IR = backsolve(chol(crossprod(W)), diag(k))
    		Btilde = crossprod(t(IR)) %*% crossprod(W, Z)
    		S = crossprod(Z - W %*% Btilde)
    		rwout = rwishart(nu + n, chol2inv(chol(V + S)))
    		B = Btilde + IR %*% matrix(rnorm(m * k), ncol = m) %*% t(rwout$CI)
    
		return(list(B = B, Sigma = rwout$IW))
	}

	rwishart = function (nu, V) {
    		m = nrow(V)
    		df = (nu + nu - m + 1) - (nu - m + 1):nu
    		if (m > 1) {
        		T = diag(sqrt(rchisq(c(rep(1, m)), df)))
        		T[lower.tri(T)] = rnorm((m * (m + 1)/2 - m))
    		}
    		else {
        		T = sqrt(rchisq(1, df))
    		}
    		U = chol(V)
    		C = t(T) %*% U
    		CI = backsolve(C, diag(m))
    
		return(list(W = crossprod(C), IW = crossprod(t(CI)), C = C, CI = CI))
	}

	getEffectsCodedParameters = function(parms, map, xcoding) {
		out = NULL
		# FOR THE NON-CODED NONE PARAMETER
		if (nrow(map) > length(xcoding)) { xcoding = c(xcoding, 1) }
		for (i in 1:nrow(map)) {
			if (xcoding[i] == 0) {
				out = cbind(out, parms[,map[i,], drop = FALSE], -1 * rowSums(parms[,map[i,], drop = FALSE]))
			}
			else {
				out = cbind(out, parms[,map[i,], drop = FALSE])
			}
		}
		return(out)
	}
	
	fsh = function() {
		if (Sys.info()[1] == "Windows") { flush.console() }
		return()
	}

	#
	# OPEN CONNECTION FOR WRITING TO LOG FILE
	#
	on.exit(sink())
	sink("RLog.txt", append = TRUE, type = "output", split = TRUE)

	#
	# PRINT TO CONSOLE AND LOG
	#
	if (restart) {
		cat("Restarting from previous ", fR, "-iteration model estmation...", sep = "", fill = TRUE)
		cat("", fill = TRUE)
	}
	cat("                    Logit Data                    ", fill = TRUE)
	cat("==================================================", fill = TRUE)
	cat("Attribute       Type         Levels", fill = TRUE)
	cat("-----------------------------------", fill = TRUE)    
	for (i in 1:natts) {
		cat(sprintf("%-12s   %-12s   %2i", paste("Attribute", i), 
			ifelse(xcoding[i] == 0, "Part Worth", "Linear"), nlevels[i]), fill = TRUE)
	} 
	cat("", fill = TRUE)
	cat(npar, " parameters to be estimated", ifelse(none, " (including 'None').", "."), sep = "", fill = TRUE)
	cat("", fill = TRUE)
	cat(nunits, "total units.", fill = TRUE)
	cat("Average of", round(avgalts, 1), "alternatives in each of", round(avgsets, 1), "sets per unit.", fill = TRUE)
	cat(ifelse(share, sum(ytab[1,]), sum(nsets)), ifelse(share, " expanded", ""), " tasks in total.", sep = "", fill = TRUE)
	cat("", fill = TRUE)
	cat("Table of choice data pooled across units:", fill = TRUE)
	cat("Choice  Count   Pct.", fill = TRUE)
	cat("--------------------", fill = TRUE)
	for (i in 1:maxalts) {
		cat(sprintf("%4i    %-6i %s", i, ytab[1,i], paste(ytab[2,i],'%', sep = "")), 
			fill = TRUE)
	}	
	cat("", fill = TRUE)
	cat("      MCMC Inference for Hierarchical Logit       ", fill = TRUE)
	cat("==================================================", fill = TRUE)
	cat("Total Iterations:         ", R + fR, fill = TRUE)
	cat("Draws used in estimation: ", use, fill = TRUE)
	cat("Units:                    ", nunits, fill = TRUE)
	cat("Parameters per unit:      ", npar, fill = TRUE)
	if (share) {
		cat("Task weight:              ", wgt, fill = TRUE)
	}
	cat("Constraints ", ifelse(constrain, "", "not "), "in effect.", fill = TRUE, sep = "")
	cat("Draws ", ifelse(save, "are to be ", "not "), "saved.", fill = TRUE, sep = "")
	cat("Prior degrees of freedom: ", df, fill = TRUE)
	cat("Prior variance:           ", v, fill = TRUE)
	cat("", fill = TRUE)
	cat("MCMC Iteration Beginning...", fill = TRUE)
	fsh()

	itime = proc.time()[3]
	acceptr.t = 0

	for (rep in 1:R) {
		if (drawdelta) {
			mgout = rGibbs(oldbetas - demos %*% t(matrix(olddelta, ncol = nz)), mubar, Amu, nu, V)
			oldcomp = mgout
			olddelta = drawDelta(demos, oldbetas, oldcomp, deltabar, Ad)
		}
		else {	
			mgout = rGibbs(oldbetas, mubar, Amu, nu, V)
			oldcomp = mgout
		}
		
		if (rep == 1) {
			if (constrain) {
				oldll = getLLMnl(oldbetas.c, y, X, info)
			}
			else {
				oldll = getLLMnl(oldbetas, y, X, info)
			}
		}

		rootpi = oldcomp$rooti
		inc.root = chol(chol2inv(chol(tcrossprod(rootpi))))
		if (drawdelta) {
			betabar = t(matrix(olddelta, ncol = nz) %*% t(demos) + oldcomp$mu)
		}
		else {
			betabar = oldcomp$mu
		}		

		if (constrain) {
			metropout = mnlRwMetropOnce(y, X, oldbetas, oldll, s, inc.root, betabar, rootpi, info, constraints, oldbetas.c)
			oldbetas.c = metropout$betadraw.c
		}
		else {
			metropout = mnlRwMetropOnce(y, X, oldbetas, oldll, s, inc.root, betabar, rootpi, info)
		}
		oldbetas = metropout$betadraw
		oldll = metropout$oldll

		#
		# CALCULATE GOODNESS-OF-FIT MEASUREMENTS
		#
		if (share) {
			RLH = exp(mean(oldll))^(1 / (avgsets * wgt))
			PctCert = (mean(oldll) - log(1 / avgalts) * avgsets * wgt) / (-log(1 / avgalts) * avgsets * wgt)
		}
		else {
			RLH = exp(mean(oldll))^(1 / avgsets)
			PctCert = (mean(oldll) - log(1 / avgalts) * avgsets) / (-log(1 / avgalts) * avgsets)
		}
		AvgVar = mean(diag(crossprod(inc.root)))
		RMS = sqrt(mean(oldbetas ^ 2))

		if (rep == 1) {
			RLH.a = RLH
			PctCert.a = PctCert
			AvgVar.a = AvgVar
			RMS.a = RMS
		}
		else {
			RLH.a = 0.99 * RLH.a + 0.01 * RLH
			PctCert.a = 0.99 * PctCert.a + 0.01 * PctCert
			AvgVar.a = 0.99 * AvgVar.a + 0.01 * AvgVar
			RMS.a = 0.99 * RMS.a + 0.01 * RMS
		}

		#
		# ADJUST SCALING PARAMETER TO TRY AND KEEP ACCEPTANCE RATE ~30%
		#
		acceptr = nunits - metropout$stay
		acceptr.t = acceptr.t + acceptr
		if (adjust.s) {
			if (acceptr / nunits < 0.3) { s = s * 0.9 }
			else if (acceptr / nunits > 0.3) { s = s * 1.1 }
		}
		acceptr = 0
		 
		#
		# PREPARE MU VALUES FOR PLOTTING
		#
		mutemp = matrix(oldcomp$mu, nrow = 1)
		muplot[rep,] = getEffectsCodedParameters(mutemp, effectsmap, xcoding)

		newplot = FALSE
		if (rep == 1) {
			yl = c(-1, 1)
			matplot(0, 0, type = "l", col = 1:8, lty = 1, ylim = yl, xlim = c(1, R + fR), xlab = "Rep", ylab = "Mu")
		}
		if (rep %% 100 == 0) {
			#
			# PLOT MU
			#
			if (max(muplot[(rep - 99):rep,]) >= yl[2]) {
				yl[2] = ceiling(max(muplot[(rep - 99):rep,]))
				newplot = TRUE
			}
			if (min(muplot[(rep - 99):rep,]) <= yl[1]) {
				yl[1] = floor(min(muplot[(rep - 99):rep,]))
				newplot = TRUE
			}
			if (newplot) {
				matplot(xplot[1:rep], muplot[1:rep,], type = "l", col = 1:8, lty = 1, ylim = yl, xlim = c(1, R + fR), xlab = "Rep", ylab = "Mu")
			}
			else {
				matplot(xplot[(rep - 100):rep], muplot[(rep - 100):rep,], type = "l", col = 1:8, lty = 1, add = TRUE)
			}

			#
			# UPDATE PROGRESS
			#
			ctime = proc.time()[3]
			timetoend = ((ctime - itime) / rep) * (R + 1 - rep)
			if (rep == 100) {
				cat("Iteration ", "Acceptance  ", "RLH    ", "Pct. Cert.  ", "Avg. Var.  ", "RMS    ", "Time to End", fill = TRUE)
			}
			cat(sprintf("%9.0i  %-5.3f        %-5.3f   %-5.3f        %-5.2f       %-5.2f   %-6s",
				rep + fR, acceptr.t / (100 * nunits), RLH.a, PctCert.a, AvgVar.a, RMS.a, 
				paste(timetoend %/% 60, ifelse(round(timetoend %% 60) > 9, ":", ":0"), round(timetoend %% 60), 
				sep = "")), fill = TRUE)
			fsh()
			acceptr.t = 0
		}
		
		if (rep > R - use) {
			if (save) {
				mkeep = (rep - (R - use)) / keep
				if (rep %% keep == 0) {
					betadraw[,,mkeep] = oldbetas
					if (constrain) { betadraw.c[,,mkeep] = oldbetas.c }
					if (drawdelta) { deltadraw[mkeep,] = olddelta }
					loglike[mkeep] = sum(oldll)
					compdraw[[mkeep]] = oldcomp
				}
			}
			if (constrain) {
				betaout = betaout + oldbetas.c
			}
			else {
				betaout = betaout + oldbetas
			}
		}		
	}
	
	ctime = proc.time()[3]
	cat("", fill = TRUE)
	cat("Total Time Elapsed: ", (ctime - itime) %/% 60, ifelse(round((ctime - itime) %% 60) > 9, ":", ":0"), 
		round((ctime - itime) %% 60), fill = TRUE, sep = "")
	cat("", fill = TRUE)

	#
	# WRITE RESTART FILE
	#	
	write.table(cbind(R + fR, nunits, npar, s, as.numeric(constrain), as.numeric(drawdelta)), "restart.txt", sep = " ", row.names = FALSE, col.names = FALSE)
	write.table(oldbetas, "restart.txt", sep = " ", row.names = FALSE, col.names = FALSE, append = TRUE)
	if (constrain) {
		write.table(oldbetas.c, "restart.txt", sep = " ", row.names = FALSE, col.names = FALSE, append = TRUE)
	}
	if (drawdelta) {
		write.table(t(olddelta), "restart.txt", sep = " ", row.names = FALSE, col.names = FALSE, append = TRUE)
	}
		
	#
	#
	# WRITE RESPONDENT BETAS TO CSV FILE
	#
	betaout = betaout / use
	betawrite = cbind(matrix(ID, ncol = 1), getEffectsCodedParameters(betaout, effectsmap, xcoding))
	betanames = "ID"
	for (i in 1:natts) {
		for (j in 1:nlevels[i]) {
			betanames = c(betanames, paste("A", i, "B", j, sep = ""))
		}
	}
	if (none) { betanames = c(betanames, "NONE") }
	cat("Writing estimated unit-level betas to Rbetas.csv in the working directory", fill = TRUE)
      cat("", fill=TRUE)
      write.table(betawrite, file = "RBetas.csv", sep = ",", col.names = betanames, row.names = FALSE, qmethod = "double") 			
	
	if (save) { 
		switch(1 + 1 * constrain + 2 * drawdelta,
			return(list(betadraw = betadraw, compdraw = compdraw, loglike = loglike)),
			return(list(betadraw = betadraw, betadraw.c = betadraw.c, compdraw = compdraw, loglike = loglike)),
			return(list(betadraw = betadraw, deltadraw = deltadraw, compdraw = compdraw, loglike = loglike)),
			return(list(betadraw = betadraw, betadraw.c = betadraw.c, deltadraw = deltadraw, compdraw = compdraw, loglike = loglike))
		)
	}
	else { return(NULL) }		 
}

