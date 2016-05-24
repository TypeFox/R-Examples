ictest <-
function(L, R, S, group, model = "GPH", type = "permutation", fuzz = 1e-12,
	output.scores = FALSE)
{
######################################################
# References: Fay, 1996, Biometrics, 811-822
# Fay, 1999, Statistics in Medicine, 273-285
# see help file for more info
######################################################
##############################################################
# for notation see Fay, 1999,  "Comparing Several Score Tests
#                 for Interval Censored Data"
#              Statistics in Medicine 18, 273-285
#############################################################
	n <- length(L)
	if(n != length(R))
		stop("length of the two interval vectors must be the\nsame")
	ss <- sort(unique(c(L, R, 0, Inf)))
	if(ss[1] < 0)
		stop("Use positive values for L and R, \n  try a  transformation"
			)
	m <- length(ss) - 2
	if(length(S) == m) {
		print("Assumed S vector missing values \n  for 0 and Inf")
		S <- c(1, S, 0)
	}
	else if(length(S) == (m + 1)) {
#  print("Assumed S vector missing values for Inf")
## no need to print this because the default from
## icfit does not output surv values for Inf
		S <- c(S, 0)
	}
	else if(length(S) != (m + 2))
		stop("The survival vector \nshould have a value for each \nunique value of L and R\n(except maybe for values at 0 and Inf)"
			)
	if(any(L == R)) {
		exacts <- sort(unique(R[R == L]))
		if(exacts[1] == 0)
			stop("L[i]==R[i]=0 for some i")
		for(j in 1:length(exacts))
			L[R == L & L == exacts[j]] <- ss[(1:(m + 2))[ss ==
				exacts[j]] - 1]
	}
	cc <- rep(NA, n)
	gg <- rep(0, n)
	for(i in 1:n) {
		gg[i] <- S[ss == L[i]] - S[ss == R[i]]
	}
	model.int <- charmatch(model, c("GPH", "Sun", "PO"))
	### allow some other values of model
	if(is.na(model.int)) {
		model.int <- charmatch(model, c("GPH", "SUN", "PO"))
		if(is.na(model.int))
			model.int <- charmatch(model, c("GPH", "SUN", "P0"))
		if(is.na(model.int))
			model.int <- charmatch(model, c("gph", "sun", "po"))
		if(is.na(model.int))
			stop("model is a character vector, either  GPH, Sun, or PO"
				)
	}
### define derivative function, dpdb without the x
### value
	dpdb.nox <- switch(model.int,
		function(P)
		{
			if(P > 0)
				P * log(P)
			else 0
		}
		,
		function(P)
		{
			np <- length(P)
			if(P[np] == 0)
				0
			else  - P[np] * sum((c(1, P[ - np]) - P)/c(1, P[ - np])
				  )
		}
		,
		function(P)
 - P * (1 - P))	## calculate the scores
	if(model.int == 2)
		for(i in 1:n)
			cc[i] <- (dpdb.nox(S[ss <= L[i]]) - dpdb.nox(S[ss <= R[
				i]]))
	else for(i in 1:n)
			cc[i] <- (dpdb.nox(S[ss == L[i]]) - dpdb.nox(S[ss == R[
				i]]))
	cc <- cc/gg
	cc
	x <- as.vector(group)
	if(length(x) != n)
		stop("group must be a vector with same length as L")
	if(type == "permutation") {
		if(is.numeric(x)) {
			ccbar <- mean(cc)
			xbar <- mean(x)
			V <- (1/(n - 1)) * sum((cc - ccbar)^2) * sum((x - xbar)^
				2)
			U <- sum(x * cc)
			chisq.value <- ((sum(x * cc) - n * ccbar * xbar)^2)/V
			pvalue <- 1 - pchisq(chisq.value, 1)
			df <- 1
			if(length(unique(x)) == 2)
				test <- "2-sample"
			else test <- "correlation"
		}
### now do the  k-sample test
### (see e.g., fay and shih, 1998, JASA, p. 389, eq 4)
		if(is.character(x)) {
			ux <- unique(x)
			nx <- length(ux)
			chisq.value <- 0
			for(i in 1:nx) {
				ccbari <- mean(cc[x == ux[i]])
				nni <- length(cc[x == ux[i]])
				chisq.value <- chisq.value + nni * ccbari^2
			}
			chisq.value <- ((n - 1)/sum(cc^2)) * chisq.value
			pvalue <- 1 - pchisq(chisq.value, nx - 1)
			test <- paste(nx, "-sample", sep = "")
	### create U so we can know which treatment is better
			U <- rep(NA, nx)
			for(ii in 1:nx)
				U[ii] <- sum(cc[x == ux[ii]])
			Unames <- rep("Group=", nx)
			Unames <- paste(Unames, ux)
			names(U) <- Unames
			df <- nx - 1
		}
		out <- list(scores = cc, U = U, chisq.value = chisq.value, df
			 = df, pvalue = pvalue, test = test, model = model,
			type = "permutation")
		if(output.scores == FALSE)
			out <- out[-1]
		return(out)
	}
### end of permutation part, the rest is only for
### type="score"
#### Check to see if we need to redefine the
####left and right censor vectors
## according to Fay, Biometrics, 1996, p. 816
# first find jumps in survival distribution
	p <-  - (S - c(1, S[ - (m + 2)]))	#
## because of rounding error, some values might
## be non-zero that should be zero. Use fuzz to
## correct this
	if(any(p[-1] < fuzz)) {
### now redefine L and R and S and ss
### if S[m+1]==0 and max(R)<Inf, then we added a 0
### to S when we did not need to, so we do not need
### to print message about ad hoc adjustment
### Otherwise print it
		if(!(all(p[ - c(1, m + 2)] >= fuzz) && S[m + 1] < fuzz && max(R
			) < Inf)) print(
				"It was necessary to use the ad hoc score test\nwith the redefined interval points, see Fay, 1996"
				)	## now make adjustments
		jumps <- ss[p > fuzz]
		if(jumps[1] == 0)
			stop("cannot have a jump at time=0")
		M <- length(jumps)
		Snew <- c(1, S[p > fuzz])
		TT <- c(0, jumps)
		Lnew <- rep(NA, n)
		Rnew <- rep(NA, n)	### TT[i] here equals S_{i-1}
###  in Fay, Biometrics, 1996, p. 816
		for(i in 1:M) {
			Lnew[L >= TT[i] & L < TT[i + 1]] <- TT[i]
			Rnew[R >= TT[i] & R < TT[i + 1]] <- TT[i]
		}
		Rnew[R >= TT[M + 1]] <- TT[M + 1]
	### go back to old notation but with new values
		L <- Lnew
		R <- Rnew
		S <- Snew
		ss <- TT
		m <- length(ss) - 2
	}
##### end of redefining of L, R, S, and ss  ######
	if(is.numeric(x)) {
		x <- matrix(x, n, 1)
		if(length(unique(x)) == 2)
			test <- "2-sample"
		else test <- "correlation"
	}
	if(is.character(x)) {
		ux <- unique(x)	#	if(length(ux) == 2) {
#		xout <- matrix(0, n,
#			1)
#		xout[x == ux[1]] <-
#			1
#		x <- xout
#		test <- "2-sample"
#	}
#	else {
		nx <- length(ux)
		xout <- matrix(0, n, nx)
		for(i in 1:nx) {
			xout[x == ux[i], i] <- 1
		}
		x <- xout
		test <- paste(nx, "-sample", sep = "")	#	}
	}
	q <- dim(x)[[2]]	### Get ready to do the observed information matrix
### see Appendix II of Fay, 1999, Stat in Med, 282
### Note there are 3 typos there:
###   1. d2P/dbdgam from logistic model, in the last ratio
###      the subscripts should be ell, ell -1 not j, j-1.
###   2. d2P/dgamk.dgamell from logistic model, insert
###      left bracket between first and second indicator
###      functions and the corresponding right bracket
###      goes at the end.
###   3. d2P/dbdgam from Proportional odds model, the
###      indicator function should be I(u=ell) not
###      I(u=k).
###
### First we need the 4 functions to get derivatives
### we already defined the first one
	dpdgam <- switch(model.int,
		function(P, u, k)
		{
			if(k == u && P > 0)
				P * log(P)
			else 0
		}
		,
		function(P, u, k)
		{
			np <- length(P)
			if(k <= u && P[np] > 0)
 - P[np] * ((P[k] - P[k + 1])/P[k])
			else 0
		}
		,
		function(P, u, k)
		{
			if(k == u && P > 0)
 - P * (1 - P)
			else 0
		}
		)
	d2pdb2.nox <- switch(model.int,
		function(P)
		{
			if(P > 0)
				P * log(P) * (1 + log(P))
			else 0
		}
		,
		function(P)
		{
			np <- length(P)
			if(np > 1 && P[np] > 0)
				P[np] * (sum((P[ - np] - P[-1])/P[ - np])^2 -
				  sum(((P[ - np] - P[-1])/P[ - np]) * (P[-1]/P[ -
				  np])))
			else 0
		}
		,
		function(P)
		{
			if(P > 0)
 - P * (1 - P) * (-1 + 2 * P)
			else 0
		}
		)
	d2pdbdgam.nox <- switch(model.int,
		function(P, u, ell)
		{
			if(u == ell && P > 0)
				P * log(P) * (1 + log(P))
			else 0
		}
		,
		function(P, u, ell)
		{
			np <- length(P)
			if(ell <= u && P[np] > 0)
				P[np] * ((P[ell] - P[ell + 1])/P[ell]) * (sum((
				  P[ - np] - P[-1])/P[ - np]) - P[ell + 1]/P[
				  ell])
			else 0
		}
		,
		function(P, u, ell)
		{
			if(u == ell && P > 0)
 - P * (1 - P) * (-1 + 2 * P)
			else 0
		}
		)
	d2pdgam2 <- switch(model.int,
		function(P, u, k, ell)
		{
			if(u == k && u == ell && P > 0)
				P * log(P) * (1 + log(P))
			else 0
		}
		,
		function(P, u, k, ell)
		{
			np <- length(P)
			if(k <= u && P[np] > 0) {
				if(ell <= u && k != ell)
				  P[np] * ((P[ell] - P[ell + 1])/P[ell]) * ((P[
				    k] - P[k + 1])/P[k])
				else if(ell <= u && k == ell)
				  P[np] * ((P[ell] - P[ell + 1])/P[ell]) * ((P[
				    k] - P[k + 1])/P[k]) - P[np] * ((P[k] - P[k +
				    1])/P[k]) * (P[k + 1]/P[k])
				else 0
			}
			else 0
		}
		,
		function(P, u, k, ell)
		{
			if(u == k && u == ell && P > 0)
 - P * (1 - P) * (-1 + 2 * P)
			else 0
		}
		)
	d2L.dB2 <- matrix(0, q, q)
	d2L.dgam2 <- matrix(0, m, m)
	d2L.dBdgam <- matrix(0, q, m)
	U <- matrix(0, q, 1)
	for(i in 1:n) {
### initialize matrices to zero
		dgi.dB <- matrix(0, q, 1)
		dgi.dgam <- matrix(0, m, 1)
		d2gi.dB2 <- matrix(0, q, q)
		d2gi.dgam2 <- matrix(0, m, m)
		d2gi.dBdgam <- matrix(0, q, m)
		dgi.dB <- x[i,  ] * gg[i] * cc[i]
	### calculate u values to go into functions
		uL <- (1:(m + 2))[ss == L[i]] - 1
		uR <- (1:(m + 2))[ss == R[i]] - 1
		if(model.int == 2) {
			d2gi.dB2 <- (x[i,  ] %*% t(x[i,  ])) * (d2pdb2.nox(S[ss <=
				L[i]]) - d2pdb2.nox(S[ss <= R[i]]))
			for(k in 1:m) {
				dgi.dgam[k] <- dpdgam(S[ss <= L[i]], uL, k) -
				  dpdgam(S[ss <= R[i]], uR, k)
				d2gi.dBdgam[, k] <- t(x[i,  ]) * (d2pdbdgam.nox(
				  S[ss <= L[i]], uL, k) - d2pdbdgam.nox(S[ss <=
				  R[i]], uR, k))
				for(ell in 1:m) {
				  d2gi.dgam2[k, ell] <- d2pdgam2(S[ss <= L[i]],
				    uL, k, ell) - d2pdgam2(S[ss <= R[i]], uR, k,
				    ell)
				}
			}
		}
		else {
### model="GPH" or model="PO"
			d2gi.dB2 <- (x[i,  ] %*% t(x[i,  ])) * (d2pdb2.nox(S[ss ==
				L[i]]) - d2pdb2.nox(S[ss == R[i]]))
	### for efficiency, we do not need to loop over k in 1:m
### since most values will be zero
			loop <- c(uL, uR)
			loop <- loop[loop > 0 & loop <= m]
			if(length(loop) > 0) {
				for(k in loop) {
				  dgi.dgam[k] <- dpdgam(S[ss == L[i]], uL, k) -
				    dpdgam(S[ss == R[i]], uR, k)
				  d2gi.dBdgam[, k] <- x[i,  ] * (d2pdbdgam.nox(
				    S[ss == L[i]], uL, k) - d2pdbdgam.nox(S[ss ==
				    R[i]], uR, k))
				  d2gi.dgam2[k, k] <- d2pdgam2(S[ss == L[i]],
				    uL, k, k) - d2pdgam2(S[ss == R[i]], uR, k,
				    k)
				}
			}
		}
		d2L.dB2 <- d2L.dB2 + gg[i]^(-1) * (d2gi.dB2 - gg[i]^(-1) *
			dgi.dB %*% t(dgi.dB))
		d2L.dBdgam <- d2L.dBdgam + gg[i]^(-1) * (d2gi.dBdgam - gg[i]^(
			-1) * dgi.dB %*% t(dgi.dgam))
		d2L.dgam2 <- d2L.dgam2 + gg[i]^(-1) * (d2gi.dgam2 - gg[i]^(-1) *
			dgi.dgam %*% t(dgi.dgam))
		U <- U + x[i,  ] * cc[i]
	}
### end i= 1 to n
	V <-  - (d2L.dB2 - d2L.dBdgam %*% solve(d2L.dgam2) %*% t(d2L.dBdgam))
	if(q == 1) {
		chisq.value <- t(U) %*% solve(V) %*% U
		df <- q
	}
	else {
		svdv <- svd(V)
		index <- (svdv$d > fuzz)
		ginvV <- svdv$v[, index] %*% ((1/svdv$d[index]) * t(
			svdv$u[, index]))
		chisq.value <- t(U) %*% ginvV %*% U
		df <- q - 1
	}
	pvalue <- 1 - pchisq(chisq.value, df)
	U <- as.vector(U)
	if(is.character(as.vector(group))) {
		Unames <- rep("Group=", length(ux))
		Unames <- paste(Unames, ux)
		names(U) <- Unames
	}
#  out <- list(scores = cc, d2L.dB2 = d2L.dB2, d2L.dBdgam =
#   d2L.dBdgam, d2L.dgam2 = d2L.dgam2, U = U, V = V,
#   chisq.value = chisq.value, df = df, pvalue = pvalue,
#   test = test, model = model, type = "score")
	out <- list(scores = cc, U = U, chisq.value = chisq.value, df = df,
		pvalue = pvalue, test = test, model = model, type = "score")
	if(output.scores == FALSE)
		out <- out[-1]
	out
}

