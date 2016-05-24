".onAttach"<-
function(...)
{
wvrelease()
}

#
# Create environment for some WaveThresh functions (PsiJ, ipndacw) to store
# results for reuse. Let efficient than previous versions of WaveThresh
# but plays more nicely with the R people
#
if (!exists("WTEnv", mode="environment"))	{
	WTEnv <- new.env()
	}

"LinaMayrand3" <-
structure(list(S = structure(c(-0.0662912607362388-0.0855811337270078i, 
-0.0662912607362388+0.0855811337270078i, 0.0352266456251514+0i, 
0.332671113131273+0i, 0.110485434560398-0.0855811337270078i, 
0.110485434560398+0.0855811337270078i, -0.0854411265843329+0i, 
0.806890861720468+0i, 0.662912607362388+0.171163681667578i, 0.662912607362388-0.171163681667578i, 
-0.135010726159072+0i, 0.45987820885317+0i, 0.662912607362388+0.171163681667578i, 
0.662912607362388-0.171163681667578i, 0.45987820885317+0i, -0.135010726159072+0i, 
0.110485434560398-0.0855811337270078i, 0.110485434560398+0.0855811337270078i, 
0.806890861720468+0i, -0.0854411265843329+0i, -0.0662912607362388-0.0855811337270078i, 
-0.0662912607362388+0.0855811337270078i, 0.332671113131273+0i, 
0.0352266456251514+0i), .Dim = as.integer(c(4, 6))), W = structure(c(-0.0662912607362388+0.0855811337270078i, 
-0.0662912607362388-0.0855811337270078i, 0.332671113131273+0i, 
0.0352266456251514+0i, -0.110485434560398-0.0855811337270078i, 
-0.110485434560398+0.0855811337270078i, -0.806890861720468+0i, 
0.0854411265843329+0i, 0.662912607362388-0.171163681667578i, 
0.662912607362388+0.171163681667578i, 0.45987820885317+0i, -0.135010726159072+0i, 
-0.662912607362388+0.171163681667578i, -0.662912607362388-0.171163681667578i, 
0.135010726159072+0i, -0.45987820885317+0i, 0.110485434560398+0.0855811337270078i, 
0.110485434560398-0.0855811337270078i, -0.0854411265843329+0i, 
0.806890861720468+0i, 0.0662912607362388-0.0855811337270078i, 
0.0662912607362388+0.0855811337270078i, -0.0352266456251514+0i, 
-0.332671113131273+0i), .Dim = as.integer(c(4, 6)))), .Names = c("S", 
"W"))
"LinaMayrand4" <-
structure(list(S = structure(c(-0.0177682977370364-0.0843076215447475i, 
0.102008915752387-0.140888496674900i, 0.512949613906065+0.139761114430506i, 
0.682186908447622+0.309503739778537i, 0.261320230715269-0.0265993641984858i, 
-0.0829326081014193-0.196341989489948i, -0.0493947656694662-0.0288541287014151i, 
0.00584356522937926+0.0277267464287373i), .Dim = as.integer(c(1, 
8))), W = structure(c(-0.00584356522937926+0.0277267464287373i, 
-0.0493947656694662+0.0288541287014151i, 0.0829326081014193-0.196341989489948i, 
0.261320230715269+0.0265993641984858i, -0.682186908447622+0.309503739778537i, 
0.512949613906065-0.139761114430506i, -0.102008915752387-0.140888496674900i, 
-0.0177682977370364+0.0843076215447475i), .Dim = as.integer(c(1, 
8)))), .Names = c("S", "W"))
"LinaMayrand5" <-
structure(list(S = structure(c(0.0104924505144049+0.0205904370844365i, 
-0.0131549130788862+0.0190001547113654i, -0.0480171707489855-0.0286805385686857i, 
0.00443868969370267-0.0660029379744943i, -0.0171289081256946+0.00872852869497756i, 
-0.0407762717133288-0.0282317864304761i, -0.0457735601342806-0.0701496826501424i, 
0.109045176430938-0.153497807951817i, -0.080639704153759-0.117947473548549i, 
0.0139497502179911-0.217696442313413i, 0.342248869674118+0.0140988497709936i, 
0.423036269003173+0.0594750872271794i, 0.151379708479645-0.0942236567554891i, 
0.245969162830182-0.123232560001445i, 0.772484323772727+0.144605393302011i, 
0.642829163846022+0.350360717350611i, 0.643003234585088+0.182852164538766i, 
0.501119052917861+0.350160634132963i, 0.479618312994977+0.059046616665079i, 
0.375016379640746+0.0994046669755474i, 0.643003234585088+0.182852164538766i, 
0.501119052917861+0.350160634132963i, -0.0564771558731019-0.0836581495806555i, 
-0.0349735956831048-0.248283003884364i, 0.151379708479645-0.0942236567554891i, 
0.245969162830182-0.123232560001445i, -0.0809927427988999-0.0456676283259696i, 
-0.106064370637416-0.113222843833651i, -0.080639704153759-0.117947473548549i, 
0.0139497502179911-0.217696442313413i, 0.0450707806910314+0.0140988497709936i, 
-0.0103356606306847+0.0594750872271794i, -0.0171289081256946+0.00872852869497756i, 
-0.0407762717133288-0.0282317864304761i, 0.0142495119522009+0.00120270047413905i, 
0.0106798133845187+0.0203460275629919i, 0.0104924505144049+0.0205904370844365i, 
-0.0131549130788862+0.0190001547113654i, -0.00819760743953431-0.00489641086342034i, 
0.000541697299744814-0.00805499281231948i), .Dim = as.integer(c(4, 
10))), W = structure(c(0.0104924505144049-0.0205904370844365i, 
-0.0131549130788862-0.0190001547113654i, -0.00819760743953431+0.00489641086342034i, 
0.000541697299744814+0.00805499281231948i, 0.0171289081256946+0.00872852869497756i, 
0.0407762717133288-0.0282317864304761i, -0.0142495119522009+0.00120270047413905i, 
-0.0106798133845187+0.0203460275629919i, -0.080639704153759+0.117947473548549i, 
0.0139497502179911+0.217696442313413i, 0.0450707806910314-0.0140988497709936i, 
-0.0103356606306847-0.0594750872271794i, -0.151379708479645-0.0942236567554891i, 
-0.245969162830182-0.123232560001445i, 0.0809927427988999-0.0456676283259696i, 
0.106064370637416-0.113222843833651i, 0.643003234585088-0.182852164538766i, 
0.501119052917861-0.350160634132963i, -0.0564771558731019+0.0836581495806555i, 
-0.0349735956831048+0.248283003884364i, -0.643003234585088+0.182852164538766i, 
-0.501119052917861+0.350160634132963i, -0.479618312994977+0.059046616665079i, 
-0.375016379640746+0.0994046669755474i, 0.151379708479645+0.0942236567554891i, 
0.245969162830182+0.123232560001445i, 0.772484323772727-0.144605393302011i, 
0.642829163846022-0.350360717350611i, 0.080639704153759-0.117947473548549i, 
-0.0139497502179911-0.217696442313413i, -0.342248869674118+0.0140988497709936i, 
-0.423036269003173+0.0594750872271794i, -0.0171289081256946-0.00872852869497756i, 
-0.0407762717133288+0.0282317864304761i, -0.0457735601342806+0.0701496826501424i, 
0.109045176430938+0.153497807951817i, -0.0104924505144049+0.0205904370844365i, 
0.0131549130788862+0.0190001547113654i, 0.0480171707489855-0.0286805385686857i, 
-0.00443868969370267-0.0660029379744943i), .Dim = as.integer(c(4, 
10)))), .Names = c("S", "W"))
"comp.theta" <-
function(djk, Sigma.inv)
{
	#
	# Takes in the complex wavelet coefficient d_{j,k} and the inverse 
	# of the covariance matrix Sigma.  Returns the scalar statistic
	# theta_{j,k}; this is \chi^2_2 if the coefficient contains 
	# only noise.
	#
	if(!is.complex(djk)) stop(
			"comp.theta should only be used on complex wavelet coefficients."
			)
	tmp <- cbind(Re(djk), Im(djk))
	tmp <- diag(tmp %*% Sigma.inv %*% t(tmp))
	return(tmp)
}
"cthr.negloglik" <-
function(parvec, dstarvec, Sigma, Sigma.inv, twopirtdetS, code)
{
	#
	# Compute -log likelihood of sample dstar from 
	# mixture of bivariate normal distributions.
	#
	# Each row of dstarvec should contain one coefficient.
	#
	if(code == "C") {
		SigVec <- c(Sigma[1, 1], Sigma[1, 2], Sigma[2, 2])
		di <- dstarvec[, 2]
		dr <- dstarvec[, 1]
		pnd <- length(di)
		pans <- 0
		Cout <- .C("Ccthrnegloglik",
			parvec = as.double(parvec),
			SigVec = as.double(SigVec),
			di = as.double(di),
			dr = as.double(dr),
			pnd = as.integer(pnd),
			pans = as.double(pans), PACKAGE = "wavethresh")
		return(Cout$pans)
	}
	else {
		p <- parvec[1]
		tmp <- parvec[3] * sqrt(parvec[2] * parvec[4])
		V <- matrix(c(parvec[2], tmp, tmp, parvec[4]), byrow = TRUE, ncol
			 = 2)
		VpS <- V + Sigma
		detVpS <- VpS[1, 1] * VpS[2, 2] - VpS[1, 2] * VpS[2, 1]
		VpS.inv <- matrix(c(VpS[2, 2],  - VpS[1, 2],  - VpS[2, 1],
			VpS[1, 1]), ncol = 2, byrow = TRUE)/detVpS
		twopirtdetVpS <- 2 * pi * sqrt(detVpS)
		tmp <- apply(dstarvec, 1, cthreb.mixden, p = p, twopirtdetS = 
			twopirtdetS, twopirtdetVpS = twopirtdetVpS, Sigma.inv
			 = Sigma.inv, VpS.inv = VpS.inv)
		return( - sum(log(tmp)))
	}
}
"cthreb.mixden" <-
function(dstar, p, twopirtdetS, twopirtdetVpS, Sigma.inv, VpS.inv)
{
	#
	# Compute density fn. of dstar from normal mixture
	#
	den1 <- exp(-0.5 * t(dstar) %*% VpS.inv %*% dstar)/twopirtdetVpS
	den2 <- exp(-0.5 * t(dstar) %*% Sigma.inv %*% dstar)/twopirtdetS
	return(p * den1 + (1 - p) * den2)
}
"cthreb.odds" <-
function(coefs, p, V, Sig, code = "NAG")
{
	#
	# Takes in coefs from a given level with EB-chosen prior parameters
	# p and V and DWT covariance matrix Sig.
	#
	# Returns posterior weights of coefficients being non-zero.
	#
	if(code == "C" || code == "NAG") {
		dr <- coefs[, 1]
		di <- coefs[, 2]
		nd <- length(dr)
		SigVec <- c(Sig[1, 1], Sig[1, 2], Sig[2, 2])
		VVec <- c(V[1, 1], V[1, 2], V[2, 2])
		pp <- p
		ans <- rep(0, nd)
		odds <- rep(0, nd)
		Cout <- .C("Ccthrcalcodds",
			pnd = as.integer(nd),
			dr = as.double(dr),
			di = as.double(di),
			VVec = as.double(VVec),
			SigVec = as.double(SigVec),
			pp = as.double(p),
			ans = as.double(ans),
			odds = as.double(odds),PACKAGE = "wavethresh")
		ptilde <- Cout$ans
	}
	else {
		VpS <- V + Sig
		detS <- Sig[1, 1] * Sig[2, 2] - Sig[1, 2]^2
		detVpS <- VpS[1, 1] * VpS[2, 2] - VpS[1, 2]^2
		mat <- solve(Sig) - solve(V + Sig)
		odds <- apply(coefs, 1, odds.matrix.mult, mat = mat)
		# Take care of excessively huge odds giving NAs in exp(odds/2)
		odds[odds > 1400] <- 1400
		odds <- p/(1 - p) * sqrt(detS/detVpS) * exp(odds/2)
		ptilde <- odds/(1 + odds)
	}
	if(any(is.na(ptilde))) {
		print("NAs in ptilde; setting those values to one")
		ptilde[is.na(ptilde)] <- 1
	}
	return(ptilde)
}
"cthreb.thresh" <-
function(coefs, ptilde, V, Sig, rule, code)
{
	#
	# Takes in coefs from a given level with EB-chosen 
	# prior covariance matrix V, posterior weights ptilde 
	# and DWT covariance matrix Sig.
	#
	# Returns thresholded coefficients; how the thresholding is
	# done depends on rule:
	#	rule == "hard": ptilde < 1/2 -> zero, otherwise
	#			keep unchanged (kill or keep).
	#	rule == "soft": ptilde < 1/2 -> zero, otherwise
	#			use posterior mean (kill or shrink).
	#	rule == "mean": use posterior mean (no zeros).
	#
	if(rule == "hard") {
		coefs[ptilde <= 0.5,  ] <- 0
		return(coefs)
	}
	else if(code == "C" || code == "NAG") {
		nd <- length(coefs[, 1])
		dr <- coefs[, 1]
		di <- coefs[, 2]
		ansr <- rep(0, nd)
		ansi <- rep(0, nd)
		VVec <- c(V[1, 1], V[1, 2], V[2, 2])
		SigVec <- c(Sig[1, 1], Sig[1, 2], Sig[2, 2])
		Cout <- .C("Cpostmean",
			pnd = as.integer(nd),
			dr = as.double(dr),
			di = as.double(di),
			VVec = as.double(VVec),
			SigVec = as.double(SigVec),
			ptilde = as.double(ptilde),
			ansr = as.double(ansr),
			ansi = as.double(ansi),PACKAGE = "wavethresh")
		coefs <- cbind(Cout$ansr, Cout$ansi)
	}
	else {
		stop("Unknown code or rule")
	}
	if(rule == "mean")
		return(coefs)
	coefs[ptilde <= 0.5,  ] <- 0
	return(coefs)
}
"cthresh" <-
function(data, j0 = 3, dwwt = NULL, dev = madmad, rule = "hard", filter.number
	 = 3.1, family = "LinaMayrand", plotfn = FALSE, TI = FALSE,
	details = FALSE, policy = "mws", code = "NAG", tol = 0.01)
{
	#
	# Limited parameter checking
	#
	n <- length(data)
	nlevels <- IsPowerOfTwo(n)
	if(is.na(nlevels))
		stop("Data should be of length a power of two.")
	if((rule != "hard") & (rule != "soft") & (rule != "mean")) {
		warning(paste("Unknown rule", rule, "so hard thresholding used"
			))
		rule <- "hard"
	}
	if((policy != "mws") & (policy != "ebayes")) {
		warning(paste("Unknown policy", policy, 
			"so using multiwavelet style thresholding"))
		policy <- "mws"
	}
	#
	# If 5 vanishing moments is called for, average over all 
	# Lina-Mayrand wavelets with 5 vanishing moments by recursively
	# calling cthresh; if filter.number=0 use all LimaMayrand wavelets
	#
	if(filter.number == 3 & ((family == "LinaMayrand") || (family = 
		"Lawton"))) {
		filter.number <- 3.1
		family <- "LinaMayrand"
	}
	else if(filter.number == 4 & family == "LinaMayrand")
		filter.number <- 4.1
	else if((filter.number == 5) & (family == "LinaMayrand")) {
		est1 <- cthresh(data, j0 = j0, dev = dev, rule = rule, 
			filter.number = 5.1, TI = TI, policy = 
			policy, details = FALSE, plotfn = FALSE, code = code, tol = tol
			)
		est2 <- cthresh(data, j0 = j0, dev = dev, rule = rule, 
			filter.number = 5.2, TI = TI, policy = 
			policy, details = FALSE, plotfn = FALSE, code = code, tol = tol
			)
		est3 <- cthresh(data, j0 = j0, dev = dev, rule = rule, 
			filter.number = 5.3, TI = TI, policy = 
			policy, details = FALSE, plotfn = FALSE, code = code, tol = tol
			)
		est4 <- cthresh(data, j0 = j0, dev = dev, rule = rule, 
			filter.number = 5.4, TI = TI, policy = 
			policy, details = FALSE, plotfn = FALSE, code = code, tol = tol
			)
		estimate <- (est1 + est2 + est3 + est4)/4
		if(plotfn) {
			x <- (1:n)/n
			plot(x, data, ylim = range(data, Re(estimate)))
			lines(x, Re(estimate), lwd = 2, col = 2)
		}
		return(estimate)
	}
	else if((filter.number == 0) & (family == "LinaMayrand")) {
		est1 <- cthresh(data, j0 = j0, dev = dev, rule = rule, 
			filter.number = 3.1, TI = TI, policy = 
			policy, details = FALSE, plotfn = FALSE, code = code, tol = tol
			)
		est2 <- cthresh(data, j0 = j0, dev = dev, rule = rule, 
			filter.number = 4.1, TI = TI, policy = 
			policy, details = FALSE, plotfn = FALSE, code = code, tol = tol
			)
		est3 <- cthresh(data, j0 = j0, dev = dev, rule = rule, 
			filter.number = 5.1, TI = TI, policy = 
			policy, details = FALSE, plotfn = FALSE, code = code, tol = tol
			)
		est4 <- cthresh(data, j0 = j0, dev = dev, rule = rule, 
			filter.number = 5.2, TI = TI, policy = 
			policy, details = FALSE, plotfn = FALSE, code = code, tol = tol
			)
		est5 <- cthresh(data, j0 = j0, dev = dev, rule = rule, 
			filter.number = 5.3, TI = TI, policy = 
			policy, details = FALSE, plotfn = FALSE, code = code, tol = tol
			)
		est6 <- cthresh(data, j0 = j0, dev = dev, rule = rule, 
			filter.number = 5.4, TI = TI, policy = 
			policy, details = FALSE, plotfn = FALSE, code = code, tol = tol
			)
		estimate <- (est1 + est2 + est3 + est4 + est5 + est6)/6
		if(plotfn) {
			x <- (1:n)/n
			plot(x, data, ylim = range(data, Re(estimate)))
			lines(x, Re(estimate), lwd = 2, col = 2)
		}
		return(estimate)
	}
	#
	# Take required type of wavelet transform.
	#
	if(TI==TRUE) data.wd <- wst(data, filter.number = filter.number, family = 
			family) else data.wd <- wd(data, filter.number = 
			filter.number, family = family)
	#
	# Generate covariance matrices
	#
	if(is.null(dwwt)) dwwt <- make.dwwt(nlevels = nlevels, filter.number = 
			filter.number, family = family)
	sigsq <- dev(Re(accessD(data.wd, level = nlevels - 1))) + dev(Im(
		accessD(data.wd, level = nlevels - 1)))
	Sigma <- array(0, c(nlevels, 2, 2))
	Sigma[, 1:2, 1:2] <- (sigsq * Im(dwwt))/2
	Sigma[, 1, 1] <- (sigsq * (1 + Re(dwwt)))/2
	Sigma[, 2, 2] <- (sigsq * (1 - Re(dwwt)))/2
	thr.wd <- data.wd
	if(policy == "mws") {
		#
		# Do multiwavelet style universal thresholding 
		#
		if(rule == "mean") {
			warning("Can't use posterior mean with multiwavelet style thresholding.  Using soft thresholding instead"
				)
			rule <- "soft"
		}
		lambda <- 2 * log(n)
		for(j in j0:(nlevels - 1)) {
			coefs <- accessD(data.wd, level = j)
			Sigma.inv <- solve(Sigma[j + 1,  ,  ])
			thetaj <- comp.theta(coefs, Sigma.inv)
			if(rule == "hard")
				coefs[abs(thetaj) < lambda] <- 0
			else {
				k <- Re(coefs)/Im(coefs)
				thetahat <- pmax(0, thetaj - lambda)
				varr <- Sigma[j + 1, 1, 1]
				vari <- Sigma[j + 1, 2, 2]
				covar <- Sigma[j + 1, 1, 2]
				bhatsq <- (varr * vari - covar^2) * thetahat
				bhatsq <- bhatsq/(vari * k^2 - 2 * covar * k +
					varr)
				coefs <- complex(modulus = sqrt(bhatsq * (k^2 +
					1)), argument = Arg(coefs))
			}
			thr.wd <- putD(thr.wd, level = j, v = coefs)
		}
	}
	else {
		#
		# Do empirical Bayes shrinkage/thresholding.
		# Start by finding parameters:
		#
		EBpars <- find.parameters(data.wd = data.wd, dwwt = dwwt, j0 = 
			j0, code = code, tol = tol, Sigma = Sigma)
		p <- c(EBpars$pars[, 1])
		Sigma <- EBpars$Sigma
		V <- array(0, dim = c(nlevels - 1, 2, 2))
		for(i in j0:(nlevels - 1))
			V[i,  ,  ] <- matrix(EBpars$pars[i, c(2, 3, 3, 4)],
				ncol = 2)
		#
		# Do thresholding.
		#
		for(j in j0:(nlevels - 1)) {
			coefs <- accessD(data.wd, level = j)
			coefs <- cbind(Re(coefs), Im(coefs))
			ptilde <- cthreb.odds(coefs, p = p[j], V = V[j,  ,
				], Sig = Sigma[j + 1,  ,  ], code = code)
			coefs.thr <- cthreb.thresh(coefs, ptilde = ptilde,
				V = V[j,  ,  ], Sig = Sigma[j,  ,  ], rule = 
				rule, code = code)
			thr.wd <- putD(thr.wd, level = j, v = complex(real = 
				coefs.thr[, 1], imaginary = coefs.thr[, 2]))
		}
	}
	#
	# Reconstruct
	#
	if(TI) data.rec <- AvBasis(thr.wd) else data.rec <- wr(thr.wd)
	#
	# Plot data and estimate
	#
	if(plotfn) {
		x <- (1:n)/n
		plot(x, data, ylim = range(data, Re(data.rec)))
		lines(x, Re(data.rec), lwd = 2, col = 2)
	}
	#
	# Return either just the estimate or an unweildy list.
	#
	if(details == FALSE) invisible(data.rec) else if(policy == "ebayes")
		invisible(list(data = data, data.wd = data.wd, thr.wd = thr.wd,
			estimate = data.rec, Sigma = Sigma, sigsq = sigsq,
			rule = rule, EBpars = EBpars$pars, wavelet = list(
			filter.number, family)))
	else invisible(list(data = data, data.wd = data.wd, thr.wd = thr.wd,
			estimate = data.rec, Sigma = Sigma, sigsq = sigsq,
			rule = rule, wavelet = list(filter.number, family)))
}
"filter.select" <-
function(filter.number, family = "DaubLeAsymm", constant = 1)
{
	G <- NULL
	if(family == "DaubExPhase") {
		family <- "DaubExPhase"
		#
		#
		#	The following wavelet coefficients are taken from
		#	Daubechies, I (1988) Orthonormal Bases of Wavelets
		#	Communications on Pure and Applied Mathematics. Page 980
		#	or Ten Lectures on Wavelets, Daubechies, I, 1992
		#	CBMS-NSF Regional Conference Series, page 195, Table 6.1
		#
		#	Comment from that table reads:
		#		"The filter coefficients for the compactly supported wavelets
		#		with extremal phase and highest number of vanishing moments
		#		compatible with their support width".
		#
		if(filter.number == 1) {
			#
			#
			#	This is for the Haar basis. (not in Daubechies).
			#
			H <- rep(0, 2)
			H[1] <- 1/sqrt(2)
			H[2] <- H[1]
			filter.name <- c("Haar wavelet")
		}
		else if(filter.number == 2) {
			H <- rep(0, 4)
			H[1] <- 0.482962913145
			H[2] <- 0.836516303738
			H[3] <- 0.224143868042
			H[4] <- -0.129409522551
			filter.name <- c("Daub cmpct on ext. phase N=2")
		}
		else if(filter.number == 3) {
			H <- rep(0, 6)
			H[1] <- 0.33267055295
			H[2] <- 0.806891509311
			H[3] <- 0.459877502118
			H[4] <- -0.13501102001
			H[5] <- -0.085441273882
			H[6] <- 0.035226291882
			filter.name <- c("Daub cmpct on ext. phase N=3")
		}
		else if(filter.number == 4) {
			H <- rep(0, 8)
			H[1] <- 0.230377813309
			H[2] <- 0.714846570553
			H[3] <- 0.63088076793
			H[4] <- -0.027983769417
			H[5] <- -0.187034811719
			H[6] <- 0.030841381836
			H[7] <- 0.032883011667
			H[8] <- -0.010597401785
			filter.name <- c("Daub cmpct on ext. phase N=4")
		}
		else if(filter.number == 5) {
			H <- rep(0, 10)
			H[1] <- 0.160102397974
			H[2] <- 0.603829269797
			H[3] <- 0.724308528438
			H[4] <- 0.138428145901
			H[5] <- -0.242294887066
			H[6] <- -0.032244869585
			H[7] <- 0.07757149384
			H[8] <- -0.006241490213
			H[9] <- -0.012580752
			H[10] <- 0.003335725285
			filter.name <- c("Daub cmpct on ext. phase N=5")
		}
		else if(filter.number == 6) {
			H <- rep(0, 12)
			H[1] <- 0.11154074335
			H[2] <- 0.494623890398
			H[3] <- 0.751133908021
			H[4] <- 0.315250351709
			H[5] <- -0.226264693965
			H[6] <- -0.129766867567
			H[7] <- 0.097501605587
			H[8] <- 0.02752286553
			H[9] <- -0.031582039318
			H[10] <- 0.000553842201
			H[11] <- 0.004777257511
			H[12] <- -0.001077301085
			filter.name <- c("Daub cmpct on ext. phase N=6")
		}
		else if(filter.number == 7) {
			H <- rep(0, 14)
			H[1] <- 0.077852054085
			H[2] <- 0.396539319482
			H[3] <- 0.729132090846
			H[4] <- 0.469782287405
			H[5] <- -0.143906003929
			H[6] <- -0.224036184994
			H[7] <- 0.071309219267
			H[8] <- 0.080612609151
			H[9] <- -0.038029936935
			H[10] <- -0.016574541631
			H[11] <- 0.012550998556
			H[12] <- 0.000429577973
			H[13] <- -0.001801640704
			H[14] <- 0.0003537138
			filter.name <- c("Daub cmpct on ext. phase N=7")
		}
		else if(filter.number == 8) {
			H <- rep(0, 16)
			H[1] <- 0.054415842243
			H[2] <- 0.312871590914
			H[3] <- 0.675630736297
			H[4] <- 0.585354683654
			H[5] <- -0.015829105256
			H[6] <- -0.284015542962
			H[7] <- 0.000472484574
			H[8] <- 0.12874742662
			H[9] <- -0.017369301002
			H[10] <- -0.044088253931
			H[11] <- 0.013981027917
			H[12] <- 0.008746094047
			H[13] <- -0.004870352993
			H[14] <- -0.000391740373
			H[15] <- 0.000675449406
			H[16] <- -0.000117476784
			filter.name <- c("Daub cmpct on ext. phase N=8")
		}
		else if(filter.number == 9) {
			H <- rep(0, 18)
			H[1] <- 0.038077947364
			H[2] <- 0.243834674613
			H[3] <- 0.60482312369
			H[4] <- 0.657288078051
			H[5] <- 0.133197385825
			H[6] <- -0.293273783279
			H[7] <- -0.096840783223
			H[8] <- 0.148540749338
			H[9] <- 0.030725681479
			H[10] <- -0.067632829061
			H[11] <- 0.000250947115
			H[12] <- 0.022361662124
			H[13] <- -0.004723204758
			H[14] <- -0.004281503682
			H[15] <- 0.001847646883
			H[16] <- 0.000230385764
			H[17] <- -0.000251963189
			H[18] <- 3.934732e-05
			filter.name <- c("Daub cmpct on ext. phase N=9")
		}
		else if(filter.number == 10) {
			H <- rep(0, 20)
			H[1] <- 0.026670057901
			H[2] <- 0.188176800078
			H[3] <- 0.527201188932
			H[4] <- 0.688459039454
			H[5] <- 0.281172343661
			H[6] <- -0.249846424327
			H[7] <- -0.195946274377
			H[8] <- 0.127369340336
			H[9] <- 0.093057364604
			H[10] <- -0.071394147166
			H[11] <- -0.029457536822
			H[12] <- 0.033212674059
			H[13] <- 0.003606553567
			H[14] <- -0.010733175483
			H[15] <- 0.001395351747
			H[16] <- 0.001992405295
			H[17] <- -0.000685856695
			H[18] <- -0.000116466855
			H[19] <- 9.358867e-05
			H[20] <- -1.3264203e-05
			filter.name <- c("Daub cmpct on ext. phase N=10")
		}
		else {
			stop("Unknown filter number for Daubechies wavelets with extremal phase and highest number of vanishing moments..."
				)
		}
	}
	else if(family == "DaubLeAsymm") {
		family <- "DaubLeAsymm"
		#
		#
		#       The following wavelet coefficients are taken from
		#       Ten Lectures on Wavelets, Daubechies, I, 1992
		#       CBMS-NSF Regional Conference Series, page 198, Table 6.3
		#
		#       Comment from that table reads:
		# 		"The low pass filter coefficients for the "least-asymmetric"
		#		compactly supported wavelets with maximum number of
		#		vanishing moments, for N = 4 to 10
		#
		if(filter.number == 4) {
			H <- rep(0, 8)
			H[1] <- -0.107148901418
			H[2] <- -0.041910965125
			H[3] <- 0.703739068656
			H[4] <- 1.136658243408
			H[5] <- 0.421234534204
			H[6] <- -0.140317624179
			H[7] <- -0.017824701442
			H[8] <- 0.045570345896
			filter.name <- c("Daub cmpct on least asymm N=4")
			H <- H/sqrt(2)
		}
		else if(filter.number == 5) {
			H <- rep(0, 10)
			H[1] <- 0.038654795955
			H[2] <- 0.041746864422
			H[3] <- -0.055344186117
			H[4] <- 0.281990696854
			H[5] <- 1.023052966894
			H[6] <- 0.89658164838
			H[7] <- 0.023478923136
			H[8] <- -0.247951362613
			H[9] <- -0.029842499869
			H[10] <- 0.027632152958
			filter.name <- c("Daub cmpct on least asymm N=5")
			H <- H/sqrt(2)
		}
		else if(filter.number == 6) {
			H <- rep(0, 12)
			H[1] <- 0.021784700327
			H[2] <- 0.004936612372
			H[3] <- -0.166863215412
			H[4] <- -0.068323121587
			H[5] <- 0.694457972958
			H[6] <- 1.113892783926
			H[7] <- 0.477904371333
			H[8] <- -0.102724969862
			H[9] <- -0.029783751299
			H[10] <- 0.06325056266
			H[11] <- 0.002499922093
			H[12] <- -0.011031867509
			filter.name <- c("Daub cmpct on least asymm N=6")
			H <- H/sqrt(2)
		}
		else if(filter.number == 7) {
			H <- rep(0, 14)
			H[1] <- 0.003792658534
			H[2] <- -0.001481225915
			H[3] <- -0.017870431651
			H[4] <- 0.043155452582
			H[5] <- 0.096014767936
			H[6] <- -0.070078291222
			H[7] <- 0.024665659489
			H[8] <- 0.758162601964
			H[9] <- 1.085782709814
			H[10] <- 0.408183939725
			H[11] <- -0.198056706807
			H[12] <- -0.152463871896
			H[13] <- 0.005671342686
			H[14] <- 0.014521394762
			filter.name <- c("Daub cmpct on least asymm N=7")
			H <- H/sqrt(2)
		}
		else if(filter.number == 8) {
			H <- rep(0, 16)
			H[1] <- 0.002672793393
			H[2] <- -0.0004283943
			H[3] <- -0.021145686528
			H[4] <- 0.005386388754
			H[5] <- 0.069490465911
			H[6] <- -0.038493521263
			H[7] <- -0.073462508761
			H[8] <- 0.515398670374
			H[9] <- 1.099106630537
			H[10] <- 0.68074534719
			H[11] <- -0.086653615406
			H[12] <- -0.202648655286
			H[13] <- 0.010758611751
			H[14] <- 0.044823623042
			H[15] <- -0.000766690896
			H[16] <- -0.004783458512
			filter.name <- c("Daub cmpct on least asymm N=8")
			H <- H/sqrt(2)
		}
		else if(filter.number == 9) {
			H <- rep(0, 18)
			H[1] <- 0.001512487309
			H[2] <- -0.000669141509
			H[3] <- -0.014515578553
			H[4] <- 0.012528896242
			H[5] <- 0.087791251554
			H[6] <- -0.02578644593
			H[7] <- -0.270893783503
			H[8] <- 0.049882830959
			H[9] <- 0.873048407349
			H[10] <- 1.015259790832
			H[11] <- 0.337658923602
			H[12] <- -0.077172161097
			H[13] <- 0.000825140929
			H[14] <- 0.042744433602
			H[15] <- -0.016303351226
			H[16] <- -0.018769396836
			H[17] <- 0.000876502539
			H[18] <- 0.001981193736
			filter.name <- c("Daub cmpct on least asymm N=9")
			H <- H/sqrt(2)
		}
		else if(filter.number == 10) {
			H <- rep(0, 20)
			H[1] <- 0.001089170447
			H[2] <- 0.000135245020
			H[3] <- -0.01222064263
			H[4] <- -0.002072363923
			H[5] <- 0.064950924579
			H[6] <- 0.016418869426
			H[7] <- -0.225558972234
			H[8] <- -0.100240215031
			H[9] <- 0.667071338154
			H[10] <- 1.0882515305
			H[11] <- 0.542813011213
			H[12] <- -0.050256540092
			H[13] <- -0.045240772218
			H[14] <- 0.07070356755
			H[15] <- 0.008152816799
			H[16] <- -0.028786231926
			H[17] <- -0.001137535314
			H[18] <- 0.006495728375
			H[19] <- 8.0661204e-05
			H[20] <- -0.000649589896
			filter.name <- c("Daub cmpct on least asymm N=10")
			H <- H/sqrt(2)
		}
		else {
			stop("Unknown filter number for Daubechies wavelets with\n least asymmetry and highest number of vanishing moments..."
				)
		}
	}
	else if(family == "MagKing") {
		family <- "MagKing"
		if(filter.number == 4) {
			H <- c(1-1i, 4-1i, 4+1i, 1+1i)/10
			G <- c(-1-2i, 5+2i, -5+2i, 1-2i)/14
			filter.name <- c("MagareyKingsbury Wavelet 4-tap")
		}
		else stop("Only have 4-tap filter at present")
	}
	else if(family == "Nason") {
		family <- "Nason"
		if(filter.number == 3) {
			H <- c(-0.066291+0.085581i,
				0.110485+0.085558i, 
				0.662912-0.171163i, 
				0.662912-0.171163i, 
				0.110485+0.085558i, 
				-0.066291+0.085581i)
			G <- c(-0.066291+0.085581i,
				-0.110485-0.085558i, 
				0.662912-0.171163i, 
				-0.662912+0.171163i
				, 0.110485+0.085558i, 
				0.066291-0.085581i)
			filter.name <- c("Nason Complex Wavelet 6-tap")
		}
		else stop("Only have 6-tap filter at present")
	}
	else if(family == "Lawton") {
		family <- "Lawton"
		if(filter.number == 3) {
			H <- c(-0.066291+0.085581i,
				0.110485+0.085558i,
				0.662912-0.171163i,
				0.662912-0.171163i,
				0.110485+0.085558,
				-0.066291+0.085581i)
			G <- c(-0.066291-0.085581i,
				-0.110485+0.085558i, 
				0.662912+0.171163i, 
				-0.662912-0.171163i
				, 0.110485-0.085558i, 
				0.066291+0.085581i)
			filter.name <- c("Lawton Complex Wavelet 6-tap")
		}
		else stop("Only have 6-tap filter at present")
	}
	else if(family == "LittlewoodPaley") {
		family <- "LittlewoodPaley"
		#
		#
		#		Define the function that computes the coefficients
		#
		hn <- function(n)
		{
			if(n == 0)
				return(1)
			else {
				pin2 <- (pi * 1:n)/2
				pin2 <- (sin(pin2)/pin2)
				return(c(rev(pin2), 1, pin2))
			}
		}
		# Next line changed in 4.6.4: added division by sqrt(2)
		H <- hn(filter.number)/sqrt(2)
		filter.name <- paste("Littlewood-Paley, N=", filter.number)
	}
	else if(family == "Yates") {
		if(filter.number != 1)
			stop("Only filter number 1 exists for Yates wavelet")
		family <- "Yates"
		H <- c(-1, 1)/sqrt(2)
		filter.name <- "Yates"
	}
	else if(family == "LinaMayrand") {
		origfn <- filter.number
		nsolution <- as.character(filter.number)
		dotpos <- regexpr("\\.", nsolution)
		leftint <- substring(nsolution, first = 1, last = dotpos - 1)
		rightint <- substring(nsolution, first = dotpos + 1, last = 
			nchar(nsolution))
		if(nchar(nsolution) == 0)
			nsolution <- 1
		else nsolution <- as.numeric(rightint)
		filter.number <- as.numeric(leftint)
		matname <- paste(family, filter.number, sep = "")
		if(!exists(matname)) {
			stop(paste("Filter matrix \"", matname, 
				"\" does not exist", sep = ""))
		}
		else {
			fm <- get(matname)
			if(nsolution > nrow(fm$S))
				stop(paste("Solution number ", nsolution, 
					" is too big. Filter matrix ", matname,
					" only has ", nrow(fm$S), " solutions")
					)
			H <- fm$S[nsolution,  ]
			G <- fm$W[nsolution,  ]
			filter.name <- paste("Lina Mayrand, J=", filter.number,
				" (nsolution=", nsolution, ")", sep = "")
		}
		filter.number <- origfn
	}
	else {
		stop("Unknown family")
	}
	H <- H/constant
	return(list(H = H, G = G, name = filter.name, family = family, 
		filter.number = filter.number))
}
"find.parameters" <-
function(data.wd, dwwt, j0, code, tol, Sigma)
{
	#
	# Preliminaries
	#
	nlevels <- nlevelsWT(data.wd)
	pars <- matrix(0, ncol = 4, nrow = nlevels - 1)
	dimnames(pars) <- list(paste("level", 1:(nlevels - 1)), c("p", 
		"var(re)", "covar(re,im)", "var(im)"))
	lower <- c(tol, tol, tol - 1, tol)
	upper <- c(1 - tol, 1000, 1 - tol, 1000)
	#
	# Calculate the covariance matrix of white noise put
	# through the DWT:
	#
	detSigma <- rep(0, nlevels)
	Sigma.inv <- array(0, c(nlevels, 2, 2))
	for(i in 1:nlevels) {
		detSigma[i] <- Sigma[i, 1, 1] * Sigma[i, 2, 2] - Sigma[i, 1,
			2]^2
		Sigma.inv[i,  ,  ] <- solve(Sigma[i,  ,  ])
	}
	#
	# Now search at each level in turn.
	#
	for(j in j0:(nlevels - 1)) {
		#
		# Get a starting point for the 
		# search over p_j and V_j 
		#
		coefs <- accessD(data.wd, level = j)
		re <- Re(coefs)
		im <- Im(coefs)
		start <- c(min(1 - 10 * tol, 0.5^(j - j0)), var(re), cor(re,
			im), var(im))
		#
		# Find the MML parameter values
		#
		coefs <- accessD(data.wd, level = j)
		dstarvec <- cbind(Re(coefs), Im(coefs))
		if(code == "NAG") {
			write(c(Sigma[j + 1, 1, 1], Sigma[j + 1, 1, 2], Sigma[
				j + 1, 2, 2]), file = "cthresh.maxloglik.data")
			write(length(re), file = "cthresh.maxloglik.data",
				append = TRUE)
			write(t(cbind(re, im)), file = "cthresh.maxloglik.data",
				append = TRUE, ncolumns = 2)
			write(start, file = "cthresh.maxloglik.start")
			write(t(cbind(lower, upper)), file = 
				"cthresh.maxloglik.start", append = TRUE)
			system("./cthresh.maxloglik")
			tmp <- scan(file = "cthresh.maxloglik.out", multi.line
				 = TRUE, quiet = TRUE)
			pars[j,  ] <- tmp[1:4]
			pars[j, 3] <- pars[j, 3] * sqrt(pars[j, 2] * pars[
				j, 4])
			ifail <- tmp[6]
			if(ifail > 0)
				warning(paste("At level", j, 
					"NAG routine e04jyf returned ifail",
					ifail))
			system("rm cthresh.maxloglik.out cthresh.maxloglik.data cthresh.maxloglik.start"
				)
		}
		else {
			if(exists("optim"))
			tmp <- optim(start, cthr.negloglik, method = 
				"L-BFGS-B", lower = lower,
				upper = upper, dstarvec = dstarvec, Sigma = 
				Sigma[j + 1,  ,  ], Sigma.inv = Sigma.inv[
				j + 1,  ,  ], twopirtdetS = 2 * pi * sqrt(
				detSigma[j + 1]), code = code)$par
			else
			tmp <- nlminb(start, cthr.negloglik, lower = lower,
				upper = upper, dstarvec = dstarvec, Sigma = 
				Sigma[j + 1,  ,  ], Sigma.inv = Sigma.inv[
				j + 1,  ,  ], twopirtdetS = 2 * pi * sqrt(
				detSigma[j + 1]), code = code)$parameters
			pars[j,  ] <- tmp
			pars[j, 3] <- pars[j, 3] * sqrt(pars[j, 2] * pars[
				j, 4])
		}
	}
	invisible(list(pars = pars, Sigma = Sigma))
}
"make.dwwt" <-
function(nlevels, filter.number = 3.1, family = "LinaMayrand")
{
	#
	# Given a choice of wavelet and number of 
	# resolution levels, compute the distinct 
	# elements of diag(WW^T).
	#
	zero.wd <- wd(rep(0, 2^nlevels), filter.number = filter.number, family
		 = family)
	dwwt <- rep(0, nlevels)
	tmp.wd <- putD(zero.wd, v = 1, level = 0)
	tmp <- Conj(wr(tmp.wd))
	#
	# tmp contains the row of W which gives the mother wavelet
	# coefficient.  Need Conj() as the inverse DWT corresponds to
	# Conj(W^T).  Now get the corresponding element of diag(WW^T)
	# by summing the squared elements of tmp.
	#
	# Then repeat for each resolution level.
	#
	dwwt[1] <- sum(tmp * tmp)
	for(lev in 1:(nlevels - 1)) {
		tmp.wd <- putD(zero.wd, v = c(1, rep(0, 2^lev - 1)), level = 
			lev)
		tmp <- Conj(wr(tmp.wd))
		dwwt[lev + 1] <- sum(tmp * tmp)
	}
	return(dwwt)
}
"odds.matrix.mult" <-
function(coef, mat)
{
	return(t(coef) %*% mat %*% coef)
}
"test.dataCT" <-
function(type = "ppoly", n = 512, signal = 1, rsnr = 7, plotfn = FALSE)
{
	x <- seq(0., 1., length = n + 1)[1:n]
	if(type == "ppoly") {
		y <- rep(0., n)
		xsv <- (x <= 0.5)
		y[xsv] <- -16. * x[xsv]^3. + 12. * x[xsv]^2.
		xsv <- (x > 0.5) & (x <= 0.75)
		y[xsv] <- (x[xsv] * (16. * x[xsv]^2. - 40. * x[xsv] + 28.))/
			3. - 1.5
		xsv <- x > 0.75
		y[xsv] <- (x[xsv] * (16. * x[xsv]^2. - 32. * x[xsv] + 16.))/
			3.
	}
	else if(type == "blocks") {
		t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44,
			0.65, 0.76, 0.78, 0.81)
		h <- c(4., -5., 3., -4., 5., -4.2, 2.1, 4.3, 
			-3.1, 2.1, -4.2)
		y <- rep(0., n)
		for(i in seq(1., length(h))) {
			y <- y + (h[i] * (1. + sign(x - t[i])))/2.
		}
	}
	else if(type == "bumps") {
		t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44,
			0.65, 0.76, 0.78, 0.81)
		h <- c(4., 5., 3., 4., 5., 4.2,	2.1, 4.3, 
			3.1, 5.1, 4.2)
		w <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03,
			0.01, 0.01, 0.005, 0.008, 0.005)
		y <- rep(0, n)
		for(j in 1:length(t)) {
			y <- y + h[j]/(1. + abs((x - t[j])/w[j]))^4.
		}
	}
	else if(type == "heavi")
		y <- 4. * sin(4. * pi * x) - sign(x - 0.3) -
			sign(0.72 - x)
	else if(type == "doppler") {
		eps <- 0.05
		y <- sqrt(x * (1. - x)) * sin((2. * pi * (1. + eps))/(x + eps))
	}
	else {
		cat(c("test.dataCT: unknown test function type", type, "\n"))
		cat(c("Terminating\n"))
		return("NoType")
	}
	y <- y/sqrt(var(y)) * signal
	ynoise <- y + rnorm(n, 0, signal/rsnr)
	if(plotfn == TRUE) {
		if(type == "ppoly")
			mlab <- "Piecewise polynomial"
		if(type == "blocks")
			mlab <- "Blocks"
		if(type == "bumps")
			mlab <- "Bumps"
		if(type == "heavi")
			mlab <- "HeaviSine"
		if(type == "doppler")
			mlab <- "Doppler"
		plot(x, y, type = "l", lwd = 2, main = mlab, ylim = range(
			c(y, ynoise)))
		lines(x, ynoise, col = 2)
		lines(x, y)
	}
	return(list(x = x, y = y, ynoise = ynoise, type = type, rsnr = rsnr))
}
"wd"<-
function(data, filter.number = 10, family = "DaubLeAsymm", type = "wavelet", bc
     = "periodic", verbose = FALSE, min.scale = 0, precond = TRUE)
{
    if(verbose == TRUE)
        cat("wd: Argument checking...")
    if(!is.atomic(data))
        stop("Data is not atomic")
    DataLength <- length(data)  #
#
# Check that we have a power of 2 data elements
#
    nlevels <- nlevelsWT(data)
    if(is.na(nlevels)) stop("Data length is not power of two")  #
#
# Check for correct type
#
    if(type != "wavelet" && type != "station")
        stop("Unknown type of wavelet decomposition")
    if(type == "station" && bc != "periodic") stop(
            "Can only do periodic boundary conditions with station"
            )   #
#
# Select the appropriate filter
#
    if(verbose == TRUE)
        cat("...done\nFilter...")
    if(bc != "interval") filter <- filter.select(filter.number = 
            filter.number, family = family) #
#
# Build the first/last database
#
    if(verbose == TRUE)
        cat("...selected\nFirst/last database...")
    fl.dbase <- first.last(LengthH = length(filter$H), DataLength = 
        DataLength, type = type, bc = bc)   #
#
#
# Check if we are doing "wavelets on the interval". If so, do it!
#
    if(bc == "interval") {
        ans <- wd.int(data = data, preferred.filter.number = 
            filter.number, min.scale = min.scale, precond = precond
            )
        fl.dbase <- first.last(LengthH = length(filter$H), DataLength
             = DataLength, type = type, bc = bc, current.scale = 
            min.scale)  #
        filter <- list(name = paste("CDV", filter.number, sep = ""), 
            family = "CDV", filter.number = filter.number)
        l <- list(transformed.vector = ans$transformed.vector, 
            current.scale = ans$current.scale, filters.used = ans$
            filters.used, preconditioned = ans$preconditioned, date
             = ans$date, nlevels = IsPowerOfTwo(length(ans$
            transformed.vector)), fl.dbase = fl.dbase, type = type, 
            bc = bc, filter = filter)
        class(l) <- "wd"
        return(l)
    }
#
#
# Save time series attribute if there is one
#
    dtsp <- tsp(data)   #
#
# Put in the data
#
    C <- rep(0, fl.dbase$ntotal)
    C[1:DataLength] <- data #
    if(verbose == TRUE)
        error <- 1
    else error <- 0
    if(verbose == TRUE) cat("built\n")  #
#
# Compute the decomposition
#
    if(verbose == TRUE)
        cat("Decomposing...\n")
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2)
    if(is.null(nbc))
        stop("Unknown boundary condition")
    ntype <- switch(type,
        wavelet = 1,
        station = 2)
    if(is.null(filter$G)) {
        wavelet.decomposition <- .C("wavedecomp",
            C = as.double(C),
            D = as.double(rep(0, fl.dbase$ntotal.d)),
            H = as.double(filter$H),
            LengthH = as.integer(length(filter$H)),
            nlevels = as.integer(nlevels),
            firstC = as.integer(fl.dbase$first.last.c[, 1]),
            lastC = as.integer(fl.dbase$first.last.c[, 2]),
            offsetC = as.integer(fl.dbase$first.last.c[, 3]),
            firstD = as.integer(fl.dbase$first.last.d[, 1]),
            lastD = as.integer(fl.dbase$first.last.d[, 2]),
            offsetD = as.integer(fl.dbase$first.last.d[, 3]),
            ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    else {
        wavelet.decomposition <- .C("comwd",
            CR = as.double(Re(C)),
            CI = as.double(Im(C)),
            LengthC = as.integer(fl.dbase$ntotal),
            DR = as.double(rep(0, fl.dbase$ntotal.d)),
            DI = as.double(rep(0, fl.dbase$ntotal.d)),
            LengthD = as.integer(fl.dbase$ntotal.d),
            HR = as.double(Re(filter$H)),
            HI = as.double( - Im(filter$H)),
            GR = as.double(Re(filter$G)),
            GI = as.double( - Im(filter$G)),
            LengthH = as.integer(length(filter$H)),
            nlevels = as.integer(nlevels),
            firstC = as.integer(fl.dbase$first.last.c[, 1]),
            lastC = as.integer(fl.dbase$first.last.c[, 2]),
            offsetC = as.integer(fl.dbase$first.last.c[, 3]),
            firstD = as.integer(fl.dbase$first.last.d[, 1]),
            lastD = as.integer(fl.dbase$first.last.d[, 2]),
            offsetD = as.integer(fl.dbase$first.last.d[, 3]),
            ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    if(verbose == TRUE)
        cat("done\n")
    error <- wavelet.decomposition$error
    if(error != 0) {
        cat("Error ", error, " occured in wavedecomp\n")
        stop("Error")
    }
    if(is.null(filter$G)) {
        l <- list(C = wavelet.decomposition$C, D = 
            wavelet.decomposition$D, nlevels = 
            nlevelsWT(wavelet.decomposition), fl.dbase = fl.dbase, 
            filter = filter, type = type, bc = bc, date = date())
    }
    else {
        l <- list(C = complex(real = wavelet.decomposition$CR, imaginary = 
            wavelet.decomposition$CI), D = complex(real = 
            wavelet.decomposition$DR, imaginary = wavelet.decomposition$DI
            ), nlevels = nlevelsWT(wavelet.decomposition), fl.dbase = 
            fl.dbase, filter = filter, type = type, bc = bc, date
             = date())
    }
    class(l) <- "wd"
    if(!is.null(dtsp))
        tsp(l) <- dtsp
    l
}
"wr.wd"<-
function(wd, start.level = 0, verbose = FALSE, bc = wd$bc, return.object = FALSE, 
    filter.number = wd$filter$filter.number, family = wd$filter$family, ...)
{
    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
    if(verbose == TRUE) cat("Argument checking...") #
#
#       Check class of wd
#
    if(verbose == TRUE)
        cat("Argument checking\n")
    ctmp <- class(wd)
    if(is.null(ctmp))
        stop("wd has no class")
    else if(ctmp != "wd")
        stop("wd is not of class wd")
    if(start.level < 0)
        stop("start.level must be nonnegative")
    if(start.level >= nlevelsWT(wd))
        stop("start.level must be less than the number of levels")
    if(is.null(wd$filter$filter.number))
        stop("NULL filter.number for wd")
    if(bc != wd$bc)
        warning("Boundary handling is different to original")
    if(wd$type == "station")
        stop("Use convert to generate wst object and then AvBasis or InvBasis"
            )
    if(wd$bc == "interval") {
        warning("All optional arguments ignored for \"wavelets on the interval\" transform"
            )
        return(wr.int(wd))
    }
    type <- wd$type
    filter <- filter.select(filter.number = filter.number, family = family)
    LengthH <- length(filter$H) #
#
#   Build the reconstruction first/last database
#
    if(verbose == TRUE)
        cat("...done\nFirst/last database...")
    r.first.last.c <- wd$fl.dbase$first.last.c[(start.level + 1):(wd$
        nlevels + 1),  ]    #
    r.first.last.d <- matrix(wd$fl.dbase$first.last.d[(start.level + 1):(wd$
        nlevels),  ], ncol = 3)
    ntotal <- r.first.last.c[1, 3] + r.first.last.c[1, 2] - r.first.last.c[
        1, 1] + 1
    names(ntotal) <- NULL
    C <- accessC(wd, level = start.level, boundary = TRUE)
    C <- c(rep(0, length = (ntotal - length(C))), C)
    Nlevels <- nlevelsWT(wd)- start.level
    error <- 0  #
#
#   Load object code
#
    if(verbose == TRUE)
        cat("...built\n")
    if(verbose == TRUE) {
        cat("Reconstruction...")
        error <- 1
    }
    ntype <- switch(type,
        wavelet = 1,
        station = 2)
    if(is.null(ntype))
        stop("Unknown type of decomposition")
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2)
    if(is.null(nbc))
        stop("Unknown boundary handling")
    if(!is.complex(wd$D)) {
        wavelet.reconstruction <- .C("waverecons",
            C = as.double(C),
            D = as.double(wd$D),
            H = as.double(filter$H),
            LengthH = as.integer(LengthH),
            nlevels = as.integer(Nlevels),
            firstC = as.integer(r.first.last.c[, 1]),
            lastC = as.integer(r.first.last.c[, 2]),
            offsetC = as.integer(r.first.last.c[, 3]),
            firstD = as.integer(r.first.last.d[, 1]),
            lastD = as.integer(r.first.last.d[, 2]),
            offsetD = as.integer(r.first.last.d[, 3]),
            ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    else {
        wavelet.reconstruction <- .C("comwr",
            CR = as.double(Re(C)),
            CI = as.double(Im(C)),
            LengthC = as.integer(length(C)),
            DR = as.double(Re(wd$D)),
            DI = as.double(Im(wd$D)),
            LengthD = as.integer(length(wd$D)),
            HR = as.double(Re(filter$H)),
            HI = as.double(Im(filter$H)),
            GR = as.double(Re(filter$G)),
            GI = as.double(Im(filter$G)),
            LengthH = as.integer(LengthH),
            nlevels = as.integer(Nlevels),
            firstC = as.integer(r.first.last.c[, 1]),
            lastC = as.integer(r.first.last.c[, 2]),
            offsetC = as.integer(r.first.last.c[, 3]),
            firstD = as.integer(r.first.last.d[, 1]),
            lastD = as.integer(r.first.last.d[, 2]),
            offsetD = as.integer(r.first.last.d[, 3]),
            ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    if(verbose == TRUE)
        cat("done\n")
    error <- wavelet.reconstruction$error
    if(error != 0) {
        cat("Error code returned from waverecons: ", error, "\n")
        stop("waverecons returned error")
    }
    fl.dbase <- wd$fl.dbase
    if(!is.complex(wd$D)) {
        l <- list(C = wavelet.reconstruction$C, D = 
            wavelet.reconstruction$D, fl.dbase = fl.dbase, nlevels
             = nlevelsWT(wd), filter = filter, type = type, bc = bc, 
            date = date())
    }
    else {
        l <- list(C = complex(real = wavelet.reconstruction$CR, imaginary = 
            wavelet.reconstruction$CI), D = complex(real = 
            wavelet.reconstruction$DR, imaginary = wavelet.reconstruction$
            DI), fl.dbase = fl.dbase, nlevels = nlevelsWT(wd), filter
             = filter, type = type, bc = bc, date = date())
    }
    class(l) <- "wd"
    if(return.object == TRUE)
        return(l)
    else return(accessC(l))
    stop("Shouldn't get here\n")
}
"wst"<-
function(data, filter.number = 10, family = "DaubLeAsymm", verbose = FALSE)
{
    if(verbose == TRUE)
        cat("Argument checking...")
    DataLength <- length(data)  #
#
# Check that we have a power of 2 data elements
#
    nlevels <- log(DataLength)/log(2)
    if(round(nlevels) != nlevels)
        stop("The length of data is not a power of 2")  #
    if(verbose == TRUE) {
        cat("There are ", nlevels, " levels\n")
    }
#
# Select the appropriate filter
#
    if(verbose == TRUE)
        cat("...done\nFilter...")
    filter <- filter.select(filter.number = filter.number, family = family)
#
#
# Compute the decomposition
#
    if(verbose == TRUE)
        cat("Decomposing...\n")
    newdata <- c(rep(0, DataLength * nlevels), data)
    Carray <- newdata
    error <- 0  #
#
#   See whether we are using complex wavelets
#
    if(is.null(filter$G)) {
        wavelet.station <- .C("wavepackst",
            Carray = as.double(Carray),
            newdata = as.double(newdata),
            DataLength = as.integer(DataLength),
            levels = as.integer(nlevels),
            H = as.double(filter$H),
            LengthH = as.integer(length(filter$H)),
            error = as.integer(error), PACKAGE  = "wavethresh")
    }
    else {
        wavelet.station <- .C("comwst",
            CaR = as.double(Re(Carray)),
            CaI = as.double(Im(Carray)),
            newdataR = as.double(Re(newdata)),
            newdataI = as.double(Im(newdata)),
            DataLength = as.integer(DataLength),
            levels = as.integer(nlevels),
            HR = as.double(Re(filter$H)),
            HI = as.double( - Im(filter$H)),
            GR = as.double(Re(filter$G)),
            GI = as.double( - Im(filter$G)),
            LengthH = as.integer(length(filter$H)),
            error = as.integer(error), PACKAGE = "wavethresh")
                }
    if(wavelet.station$error != 0)
        stop(paste("Memory error in wavepackst (or comwst). Code ", 
            wavelet.station))
    if(is.null(filter$G)) {
        wpm <- matrix(wavelet.station$newdata, ncol = DataLength, byrow
             = TRUE)
        Carray <- matrix(wavelet.station$Carray, ncol = DataLength, 
            byrow = TRUE)
    }
    else {
        newdata <- complex(real = wavelet.station$newdataR, imaginary = 
            wavelet.station$newdataI)
        Carray <- complex(real = wavelet.station$CaR, imaginary = 
            wavelet.station$CaI)
        wpm <- matrix(newdata, ncol = DataLength, byrow = TRUE)
        Carray <- matrix(Carray, ncol = DataLength, byrow = TRUE)
    }
    wp <- list(wp = wpm, Carray = Carray, nlevels = nlevels, filter = 
        filter, date = date())
    class(wp) <- "wst"
    wp
}
"AutoBasis"<-
function(wp, verbose = FALSE, zilchtol = 1e-08,entropy = Shannon.entropy)
{
    if(class(wp) != "wp") {
        stop("Can only operate on wavelet packet objects")
    }
    if(IsEarly(wp)) {
        ConvertMessage()
        stop()
    }
#
#
#   Including the original data set there are nlevels levels. Labelled
#   0,...,nlevels-1. Level nlevels-1 is the original data set.
#
    nlevels <- nlevelsWT(wp)
    for(i in 1:(nlevels - 1)) {
        NPBaseLev <- 2^(nlevels - i)
        PKLength <- 2^i
        if(verbose == TRUE) {
            cat("Base level is ", i)
            cat(" Number of packets is ", NPBaseLev, "\n")
            cat(" Packet Length is ", PKLength, "\n")
        }
        scan()
        for(j in 0:(NPBaseLev - 1)) {
            p1 <- getpacket(wp, level = (i - 1), index = 2 * j)
            p2 <- getpacket(wp, level = (i - 1), index = 2 * j + 1)
            p <- getpacket(wp, level = i, index = j)
            if(verbose == TRUE) {
                cat("Comparing: (", i, ",", j, ") with ")
                cat("(", (i - 1), ",", 2 * j, ") + (", (i - 1), 
                  ",", 2 * j + 1, ")\n")
            }
            if(is.na(p1[1]) || is.na(p2[1])) {
                if(verbose == TRUE) {
                  cat("Upper Level is not eligible for")
                  cat(" incorporation. Moving on...\n")
                }
                wp <- putpacket(wp, lev = i, index = j, packet
                   = rep(NA, length = length(p)))
            }
            else {
                e1 <- entropy(p1, zilchtol)
                e2 <- entropy(p2, zilchtol)
                e <- entropy(p, zilchtol)
                if(verbose == TRUE) {
                  cat("Entropy:", signif(e, 3), "?", signif(e1, 
                    3), "+", signif(e2, 3), "=", signif(e1 + e2,
                    3))
                }
                if(e < e1 + e2 || (is.infinite(e) && is.infinite(e1) && 
                  is.infinite(e2))) {
                  wp <- putpacket(wp, level = (i - 1), index = 
                    2 * j, packet = rep(NA, length = PKLength/2
                    ))
                  wp <- putpacket(wp, level = (i - 1), index = 
                    2 * j + 1, packet = rep(NA, length = 
                    PKLength/2))
                }
                else {
                  wp <- putpacket(wp, level = i, index = j, 
                    packet = rep(NA, length = PKLength))
                }
                if(e < e1 + e2 || (is.infinite(e) && is.infinite(e1) && 
                  is.infinite(e2)))
                  cat(" REPLACE\n")
                else cat(" KEEP\n")
            }
        }
    }
    wp
}
"AvBasis"<-
function(...)
UseMethod("AvBasis")
"AvBasis.wst"<-
function(wst, Ccode = TRUE, ...)
{
    nlevels <- nlevelsWT(wst)
    if(is.null(wst$filter$G)) {
        if(Ccode == FALSE) {
            answer <- av.basis(wst, level = nlevels - 1, ix1 = 0, 
                ix2 = 1, filter = wst$filter)
        }
        else {
            error <- 0
            answer <- rep(0, 2^nlevels)
            H <- wst$filter$H
            aobj <- .C("av_basisWRAP",
                wstR = as.double(wst$wp),
                wstC = as.double(wst$Carray),
                LengthData = as.integer(length(answer)),
                level = as.integer(nlevels - 1),
                H = as.double(H),
                LengthH = as.integer(length(H)),
                answer = as.double(answer),
                error = as.integer(error), PACKAGE = "wavethresh")
            if(aobj$error != 0)
                stop(paste("av_basisWRAP returned error code", 
                  aobj$error))
            answer <- aobj$answer
        }
    }
    else {
        error <- 0
        answerR <- answerI <- rep(0, 2^nlevels)
        H <- wst$filter$H
        G <- wst$filter$G
        aobj <- .C("comAB_WRAP",
            wstR = as.double(Re(wst$wp)),
            wstI = as.double(Im(wst$wp)),
            wstCR = as.double(Re(wst$Carray)),
            wstCI = as.double(Im(wst$Carray)),
            LengthData = as.integer(length(answerR)),
            level = as.integer(nlevels - 1),
            HR = as.double(Re(H)),
            HI = as.double(Im(H)),
            GR = as.double(Re(G)),
            GI = as.double(Im(G)),
            LengthH = as.integer(length(H)),
            answerR = as.double(answerR),
            answerI = as.double(answerI),
            error = as.integer(error), PACKAGE = "wavethresh")
        if(aobj$error != 0)
            stop(paste("av_basisWRAP returned error code", aobj$
                error))
        answer <- complex(real = aobj$answerR, imaginary = aobj$answerI)
    }
    answer
}
"AvBasis.wst2D"<-
function(wst2D, ...)
{
    filter <- wst2D$filter
    amdim <- dim(wst2D$wst2D)
    im <- matrix(0, nrow = amdim[2]/2, ncol = amdim[2]/2)
    ans <- .C("SAvBasis",
        am = as.double(wst2D$wst2D),
        d1 = as.integer(amdim[1]),
        d12 = as.integer(amdim[1] * amdim[2]),
        TheSmooth = as.double(im),
        levj = as.integer(amdim[1]),
        H = as.double(filter$H),
        LengthH = as.integer(length(filter$H)),
        error = as.integer(0), PACKAGE = "wavethresh")
    if(ans$error != 0)
        stop(paste("Error code was ", ans$error))
    matrix(ans$TheSmooth, nrow = amdim[2]/2)
}
"BAYES.THR"<-
function(data, alpha = 0.5, beta = 1, filter.number = 8, family = "DaubLeAsymm",
    bc = "periodic", dev = var, j0 = 5, plotfn = FALSE)
{
#
#------------Estimation of C1 and C2 via universal threshodling-----------------
#
    ywd <- wd(data, filter.number = filter.number, family = family, bc = bc
        )
    sigma <- sqrt(dev(accessD(ywd, level = (nlevelsWT(ywd) - 1))))
    uvt <- threshold(ywd, policy = "universal", type = "soft", dev = dev, 
        by.level = FALSE, levels = (nlevelsWT(ywd) - 1), return.threshold = TRUE)
    universal <- threshold(ywd, policy = "manual", value = uvt, type = 
        "soft", dev = dev, levels = j0:(nlevelsWT(ywd) - 1))
    nsignal <- rep(0, nlevelsWT(ywd))
    sum2 <- rep(0, nlevelsWT(ywd))
    for(j in 0:(nlevelsWT(ywd) - 1)) {
        coefthr <- accessD(universal, level = j)
        nsignal[j + 1] <- sum(abs(coefthr) > 0)
        if(nsignal[j + 1] > 0)
            sum2[j + 1] <- sum(coefthr[abs(coefthr) > 0]^2)
    }
    C <- seq(1000, 15000, 50)
    l <- rep(0, length(C))
    lev <- seq(0, nlevelsWT(ywd) - 1)
    v <- 2^( - alpha * lev)
    for(i in 1:length(C)) {
        l[i] <- 0.5 * sum(- nsignal * (log(sigma^2 + C[i] * v) + 2 * log(pnorm(( - sigma * sqrt(2 * log(2^nlevelsWT(ywd))))/
            sqrt(sigma^2 + C[i] * v)))) - sum2/2/(sigma^2 + C[i] * v))
    }
    C1 <- C[l == max(l)]
    tau2 <- C1 * v
    p <- 2 * pnorm(( - sigma * sqrt(2 * log(2^nlevelsWT(ywd))))/sqrt(sigma^2 + 
        tau2))
    if(beta == 1)
        C2 <- sum(nsignal/p)/nlevelsWT(ywd)
    else C2 <- (1 - 2^(1 - beta))/(1 - 2^((1 - beta) * nlevelsWT(ywd))) * sum(
            nsignal/p)
    pr <- pmin(1, C2 * 2^( - beta * lev))
    rat <- tau2/(sigma^2 + tau2)    #
#   
#----------------------Bayesian Thresholding------------------------------------
#
    bayesian <- ywd
    for(j in 0:(nlevelsWT(ywd)- 1)) {
        coef <- accessD(ywd, level = j)
        w <- (1 - pr[j + 1])/pr[j + 1]/sqrt((sigma^2 * rat[j + 1])/tau2[
            j + 1]) * exp(( - rat[j + 1] * coef^2)/2/sigma^2)
        z <- 0.5 * (1 + pmin(w, 1))
        median <- sign(coef) * pmax(0, rat[j + 1] * abs(coef) - sigma * 
            sqrt(rat[j + 1]) * qnorm(z))
        bayesian <- putD(bayesian, level = j, v = median)
    }
    bayesrec <- wr(bayesian)    #
#---------------Resulting plots--------------------------------------------
#
    if(plotfn == TRUE) {
        x <- seq(1, length(data))/length(data)
        par(mfrow = c(1, 2))
        plot(x, data, type = "l", ylab = "(a) Data")
        plot(x, bayesrec, type = "l", ylab = "(b) BayesThresh", ylim = 
            c(min(data), max(data)))
    }
    return(bayesrec)
}

"BMdiscr"<-
function(BP)
{
    dm <- lda(x = BP$BasisMatrix, grouping = BP$groups)   #
    BMd <- list(BP = BP, dm = dm)
}

"Best1DCols"<-
function(w2d, mincor = 0.69999999999999996)
{
    m <- w2d$m
    level <- w2d$level
    pktix <- w2d$pktix
    nbasis <- length(level)
    corvec <- rep(0, nbasis)
    for(i in 1:nbasis) {
        corvec[i] <- cor(m[, i], w2d$groups)
    }
    corvec <- abs(corvec)
    sv <- corvec > mincor
	if (sum(sv) < 2)
		stop("Not enough variables. Decrease mincor")
    m <- m[, sv]
    level <- level[sv]
    pktix <- pktix[sv]
    corvec <- corvec[sv]
    sl <- rev(sort.list(corvec))
    l <- list(nlevels = nlevelsWT(w2d), BasisMatrix = m[, sl], level = level[
        sl], pkt = pktix[sl], basiscoef = corvec[sl], groups = w2d$groups)
    class(l) <- "BP"
    l
}
"CWCV"<-
function(ynoise, ll, x = 1:length(ynoise), filter.number = 10, family = 
    "DaubLeAsymm", thresh.type = "soft", tol = 0.01, maxits=500,
	verbose = 0, plot.it
     = TRUE, interptype = "noise")
{
#
#   Switch on verbosity for function calls if necessary
#
    if(verbose == 2) CallsVerbose <- TRUE else CallsVerbose <- FALSE
    if(verbose == 1)
        cat("WaveletCV: Wavelet model building\nThinking ")
    n <- length(ynoise)
    ywd <- wd(ynoise, filter.number = filter.number, family = family, 
        verbose = CallsVerbose)
    univ.threshold <- threshold(ywd, type = thresh.type, return.threshold
         = TRUE, lev = ll:(nlevelsWT(ywd)- 1), verbose = CallsVerbose, 
        policy = "universal")[1]
    if(verbose == 1) {
        cat("Universal threshold: ", univ.threshold, "\n")
        cat("Now doing universal threshold reconstruction...")
    }
    yuvtwd <- threshold(ywd, type = thresh.type, lev = ll:(nlevelsWT(ywd)- 1),
        verbose = CallsVerbose, policy = "universal")
    if(verbose == 1)
        cat("done\nNow reconstructing...")
    yuvtwr <- wr(yuvtwd, verbose = CallsVerbose)
    if(verbose == 1)
        cat("done\nNow plotting universal thresholded\n")
    if(plot.it == TRUE) {
        oldpar <- par(mfrow = c(2, 2))
        matplot(x, cbind(ynoise, yuvtwr), type = "l", main = 
            "Universal Threshold Reconstruction", xlab = "x", col
             = c(3, 2), lty = c(3, 2))
    }
    filter <- filter.select(filter.number = filter.number, family = family)
    N <- length(ynoise)
    nlevels <- log(N)/log(2)
    ssq <- 0
    if(verbose > 0)
        error <- 1
    else error <- 0
    if(round(nlevels) != nlevels)
        stop("Datalength not power of 2")
    fl.dbase <- first.last(length(filter$H), N/2)
    C <- rep(0, fl.dbase$ntotal)
    D <- rep(0, fl.dbase$ntotal.d)
    ntt <- switch(thresh.type,
        hard = 1,
        soft = 2)
    if(is.null(ntt))
        stop("Unknown threshold type")
    interptype <- switch(interptype,
        noise = 1,
        normal = 2)
    if(is.null(interptype))
        stop("Unknown interptype")
    bc <- "periodic"
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2)
    if(is.null(nbc))
        stop("Unknown boundary conditions")
    xvthresh <- 0
    if(verbose == 1)
        cat("Now optimising cross-validated error estimate\n")
    ans <- .C("CWaveletCV",
        noisy = as.double(ynoise),
        nnoisy = as.integer(N),
        univ.threshold = as.double(univ.threshold),
        C = as.double(C),
        D = as.double(D),
        LengthD = as.integer(length(D)),
        H = as.double(filter$H),
        LengthH = as.integer(length(filter$H)),
        levels = as.integer(nlevels),
        firstC = as.integer(fl.dbase$first.last.c[, 1]),
        lastC = as.integer(fl.dbase$first.last.c[, 2]),
        offsetC = as.integer(fl.dbase$first.last.c[, 3]),
        firstD = as.integer(fl.dbase$first.last.d[, 1]),
        lastD = as.integer(fl.dbase$first.last.d[, 2]),
        offsetD = as.integer(fl.dbase$first.last.d[, 3]),
        ntt = as.integer(ntt),
        ll = as.integer(ll),
        nbc = as.integer(nbc),
        tol = as.double(tol),
	maxits = as.integer(maxits),
        xvthresh = as.double(xvthresh),
        interptype = as.integer(interptype),
        error = as.integer(error), PACKAGE = "wavethresh")

    if (ans$error == 1700)	{
		message("Algorithm not converging (yet).")
		message("Maybe increase number of maximum iterations (maxits or cvmaxits)?")
		message("Or increase tolerance (tol or cvtol) a bit?")
		message("Wanted to achieve tolerance of ", tol,
			" but have actually achieved: ", ans$tol)
		message("Check levels you are thresholding, especially if length of data set is small. E.g. if n<=16 then default levels argument probably should be changed.")
		stop(paste("Maximum number of iterations", maxits, " exceeded."))
		}
    else if(ans$error != 0) {
        cat("Error code ", ans$error, "\n")
        stop("There was an error")
    }
#
#
#   Now do the reconstuction using xvthresh
#
    xvwd <- threshold(ywd, policy = "manual", value = ans$xvthresh, type = 
        thresh.type, lev = ll:(nlevelsWT(ywd)- 1))
    xvwddof <- dof(xvwd)
    xvwr <- wr(xvwd)
    if(plot.it == TRUE)
        matplot(x, cbind(ynoise, yuvtwr, xvwr), type = "l", main = 
            "XV Threshold Reconstruction", xlab = "x", col = c(3, 2,
            1))
    fkeep <- NULL
    xkeep <- NULL
    list(x = x, ynoise = ynoise, xvwr = xvwr, yuvtwr = yuvtwr, xvthresh = 
        ans$xvthresh, uvthresh = univ.threshold, xvdof = xvwddof, uvdof
         = dof(yuvtwd), xkeep = xkeep, fkeep = fkeep)
}
"CWavDE"<-
function(x, Jmax, threshold = 0, nout = 100, primary.resolution = 1, 
    filter.number = 10, family = "DaubLeAsymm", verbose = 0, SF = NULL, WV
     = NULL)
{
    rx <- range(x)
    xout <- rep(0, nout)
    fout <- rep(0, nout)
    kmin <- 0
    kmax <- 0
    kminW <- rep(0, Jmax)
    kmaxW <- rep(0, Jmax)
    xminW <- rep(0, Jmax)
    xmaxW <- rep(0, Jmax)   #
#   Generate the scaling function and the wavelet if they're not supplied
#
    if(is.null(SF)) {
        if(verbose > 0)
            cat("Computing scaling function\n")
        SF <- draw.default(filter.number = filter.number, family = 
            family, plot.it = FALSE, scaling.function = TRUE, enhance = FALSE)
    }
    if(is.null(WV)) {
        if(verbose > 0)
            cat("Computing wavelet function\n")
        WV <- draw.default(filter.number = filter.number, family = 
            family, plot.it = FALSE, enhance = FALSE)
    }
    swv <- support(filter.number = filter.number, family = family)  #
    error <- 0
    ans <- .C("CWavDE",
        x = as.double(x),
        n = as.integer(length(x)),
        minx = as.double(rx[1]),
        maxx = as.double(rx[2]),
        Jmax = as.integer(Jmax),
        threshold = as.double(threshold),
        xout = as.double(xout),
        fout = as.double(fout),
        nout = as.integer(nout),
        primary.resolution = as.double(primary.resolution),
        SFx = as.double(SF$x),
        SFy = as.double(SF$y),
        lengthSF = as.integer(length(SF$x)),
        WVx = as.double(WV$x),
        WVy = as.double(WV$y),
        lengthWV = as.integer(length(WV$x)),
        kmin = as.integer(kmin),
        kmax = as.integer(kmax),
        kminW = as.integer(kminW),
        kmaxW = as.integer(kmaxW),
        xminW = as.double(xminW),
        xmaxW = as.double(xmaxW),
        phiLH = as.double(swv$phi.lh),
        phiRH = as.double(swv$phi.rh),
        psiLH = as.double(swv$psi.lh),
        psiRH = as.double(swv$psi.rh),
        verbose = as.integer(verbose),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(ans$error != 0)
        stop(paste("CWavDE returned error code", ans$error))
    l <- list(x = ans$xout, y = ans$fout, sfix = ans$kmin:ans$kmax, wvixmin
         = ans$kminW, wvixmax = ans$kmaxW)
    l
}
"CanUseMoreThanOneColor"<-
function()
{
#
# In the S version of this code it was possible to interrogate certain
# graphics devices to see how many colors they display.
# Most users these days will be using X11, or quartz or pdf which can
# so this routine is fixed now to return true.

return(TRUE)
}
"ConvertMessage"<-
function()
{
    cat("Your wavelet object is from an old release of wavethresh.\n")
    cat("Please apply the function convert() to your object.\n")
    cat("This will update it to the most up to date release.\n")
    cat("e.g. if the name of your wavelet object is \"fred\" then type:\n")
    cat("fred <- convert(fred)\n")
}
"Crsswav"<-
function(noisy, value = 1, filter.number = 10, family = "DaubLeAsymm", 
    thresh.type = "hard", ll = 3)
{
    filter <- filter.select(filter.number = filter.number, family = family)
    N <- length(noisy)
    nlevels <- log(N)/log(2)
    ssq <- 0
    error <- 0
    if(round(nlevels) != nlevels)
        stop("Datalength not power of 2")
    fl.dbase <- first.last(length(filter$H), N/2)
    C <- rep(0, fl.dbase$ntotal)
    D <- rep(0, fl.dbase$ntotal.d)
    ntt <- switch(thresh.type,
        hard = 1,
        soft = 2)
    if(is.null(ntt))
        stop("Unknown threshold type")
    bc <- "periodic"
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2)
    if(is.null(nbc))
        stop("Unknown boundary conditions")
    ans <- .C("Crsswav",
        noisy = as.double(noisy),
        nnoisy = as.integer(N),
        value = as.double(value),
        C = as.double(C),
        D = as.double(D),
        LengthD = as.integer(length(D)),
        H = as.double(filter$H),
        LengthH = as.integer(length(filter$H)),
        levels = as.integer(nlevels),
        firstC = as.integer(fl.dbase$first.last.c[, 1]),
        lastC = as.integer(fl.dbase$first.last.c[, 2]),
        offsetC = as.integer(fl.dbase$first.last.c[, 3]),
        firstD = as.integer(fl.dbase$first.last.d[, 1]),
        lastD = as.integer(fl.dbase$first.last.d[, 2]),
        offsetD = as.integer(fl.dbase$first.last.d[, 3]),
        ntt = as.integer(ntt),
        ll = as.integer(ll),
        nbc = as.integer(nbc),
        ssq = as.double(ssq),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(ans$error != 0) {
        cat("Error code ", ans$error, "\n")
        stop("There was an error")
    }
    cat("The answer was ", ans$ssq, "\n")
    return(list(ssq = ans$ssq, value = value, type = thresh.type, lev = ll:(
        nlevels - 1)))
}
"Cthreshold"<-
function(wd, thresh.type = "soft", value = 0, levels = 3:(nlevelsWT(wd)- 1))
{
    D <- wd$D
    Dlevels <- nlevelsWT(wd)- 1
    error <- 0
    ntt <- switch(thresh.type,
        hard = 1,
        soft = 2)
    if(is.null(ntt))
        stop("Unknown thresh.type")
    nbc <- switch(wd$bc,
        periodic = 1,
        symmetric = 2)
    if(is.null(nbc))
        stop("Unknown boundary conditions")
    ans <- .C("Cthreshold",
        D = as.double(D),
        LengthD = as.integer(wd$fl.dbase$ntotal.d),
        firstD = as.integer(wd$fl.dbase$first.last.d[, 1]),
        lastD = as.integer(wd$fl.dbase$first.last.d[, 2]),
        offsetD = as.integer(wd$fl.dbase$first.last.d[, 3]),
        Dlevels = as.integer(Dlevels),
        ntt = as.integer(ntt),
        value = as.double(value),
        levels = as.integer(levels),
        qlevels = as.integer(length(levels)),
        nbc = as.integer(nbc),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(ans$error != 0) {
        stop("Error occurred")
        cat("Error code was ", ans$error, "\n")
    }
    wd$D <- ans$D
    wd
}
"DJ.EX"<-
function(n = 1024, signal = 7, rsnr = 7, noisy = FALSE, plotfn = FALSE)
{
    x <- seq(1, n)/n    
    #--------------------Blocks---------------------------------------------------
    t <- c(0.10000000000000001, 0.13, 0.14999999999999999, 
        0.23000000000000001, 0.25, 0.40000000000000002, 0.44, 
        0.65000000000000002, 0.76000000000000001, 0.78000000000000003, 
        0.81000000000000005)
    h1 <- c(4, -5, 3, -4, 5, -4.2000000000000002, 2.1000000000000001, 
        4.2999999999999998, -3.1000000000000001, 2.1000000000000001, 
        -4.2000000000000002)
    blocks <- rep(0, n)
    for(i in seq(1, length(h1))) {
        blocks <- blocks + (h1[i] * (1 + sign(x - t[i])))/2
    }
#--------------------Bumps----------------------------------------------------
    h2 <- c(4, 5, 3, 4, 5, 4.2000000000000002, 2.1000000000000001, 
        4.2999999999999998, 3.1000000000000001, 5.0999999999999996, 
        4.2000000000000002)
    w <- c(0.0050000000000000001, 0.0050000000000000001, 
        0.0060000000000000001, 0.01, 0.01, 0.029999999999999999, 0.01, 
        0.01, 0.0050000000000000001, 0.0080000000000000002, 
        0.0050000000000000001)
    bumps <- rep(0, n)
    for(i in seq(1, length(h2))) {
        bumps <- bumps + h2[i] * pmax(0, (1 - abs((x - t[i])/w[i])))^4
    }
#-------------------HeaviSine-------------------------------------------------
    heavi <- 4 * sin(4 * pi * x) - sign(x - 0.29999999999999999) - sign(
        0.71999999999999997 - x)    
    #--------------------Doppler--------------------------------------------------
    eps <- 0.050000000000000003
    doppler <- sqrt(x * (1 - x)) * sin((2 * pi * (1 - eps))/(x + eps))  
    #------------------------Normalization----------------------------------------
    blocks <- blocks/sqrt(var(blocks)) * signal
    bumps <- bumps/sqrt(var(bumps)) * signal
    heavi <- heavi/sqrt(var(heavi)) * signal
    doppler <- doppler/sqrt(var(doppler)) * signal
    if(noisy == TRUE) {
        values <- list(blocks = blocks + rnorm(n, 0, signal/rsnr), 
            bumps = bumps + rnorm(n, 0, signal/rsnr), heavi = heavi +
            rnorm(n, 0, signal/rsnr), doppler = doppler + rnorm(n, 
            0, signal/rsnr))
    }
    else {
        values <- list(blocks = blocks, bumps = bumps, heavi = heavi, 
            doppler = doppler)
    }
    if(plotfn == TRUE) {
        par(mfrow = c(3, 2))
        plot(x, values$blocks, type = "l", ylab = "(a) Blocks")
        plot(x, values$bumps, type = "l", ylab = "(b) Bumps")
        plot(x, values$heavi, type = "l", ylab = "(c) HeaviSine")
        plot(x, values$doppler, type = "l", ylab = "(d) Doppler")
    }
    return(values)
}
"FullWaveletCV"<-
function(noisy, ll = 3, type = "soft", filter.number = 10, family = 
    "DaubLeAsymm", tol = 0.01, verbose = 0)
{
    noisywd <- wd(noisy, filter.number = filter.number, family = family)
    softuv <- threshold(noisywd, levels = ll:(nlevelsWT(noisywd)- 1), type = 
        "soft", policy = "universal", dev = madmad, return.thresh = TRUE)
    H <- filter.select(filter.number = filter.number, family = family)$H
    ntt <- switch(type,
        hard = 1,
        soft = 2)
    error <- verbose
    xvthresh <- 0
    ans <- .C("FullWaveletCV",
        noisy = as.double(noisy),
        nnoisy = as.integer(length(noisy)),
        UniversalThresh = as.double(softuv),
        H = as.double(H),
        LengthH = as.integer(length(H)),
        ntt = as.integer(ntt),
        ll = as.integer(ll),
        tol = as.double(tol),
        xvthresh = as.double(xvthresh),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(ans$error != 0) {
        cat("Error code returned was ", ans$error, "\n")
        stop("Error detected from C routine")
    }
    ans$xvthresh
}
"GenW"<-
function(n = 8, filter.number = 10, family = "DaubLeAsymm", bc = "periodic")
{
    z <- rep(0, n)
    if(bc == "periodic") {
        w <- matrix(0, nrow = n, ncol = n)
        for(i in 1:n) {
            v <- z
            v[i] <- 1
            wobj <- wd(v, filter.number = filter.number, family = 
                family, bc = bc)
            w[i, 1] <- accessC(wobj, lev = 0)
            w[i, 2:n] <- wobj$D
        }
    }
    else {
        w <- NULL
        for(i in 1:n) {
            v <- z
            v[i] <- 1
            wobj <- wd(v, filter.number = filter.number, family = 
                family, bc = bc)
            wrow <- c(accessC(wobj, lev = 0, boundary = TRUE), wobj$D)
            w <- rbind(w, wrow)
        }
    }
    w
}
"GetRSSWST"<-
function(ndata, threshold, levels, family = "DaubLeAsymm", filter.number = 10, 
    type = "soft", norm = l2norm, verbose = 0, InverseType = "average")
{
    thverb <- FALSE
    if(verbose > 1)
        thverb <- TRUE
    if(InverseType != "average" && InverseType != "minent") stop(paste(
            "Unknown InverseType: ", InverseType))  #
# Get odds and evens
#
    oddsv <- seq(from = 1, to = length(ndata), by = 2)
    evensv <- seq(from = 2, to = length(ndata), by = 2)
    odata <- ndata[oddsv]
    edata <- ndata[evensv]  #
#
# Build odd thresholded estimate, then, threshold and rebuild
#
    odataWST <- wst(odata, filter.number = filter.number, family = family)
    odataWSTt <- threshold.wst(odataWST, levels = levels, policy = "manual",
        value = threshold, verbose = thverb)
    if(InverseType == "average")
        odataWSTr <- AvBasis.wst(odataWSTt) #
    else if(InverseType == "minent") {
        odataNV <- MaNoVe(odataWSTt)
        cat("ODD Node Vector\n")
        cat("---------------\n")
        print(odataNV)
        odataWSTr <- InvBasis.wst(odataWSTt, nv = odataNV)
    }
    else stop(paste("Unknown InverseType: ", InverseType))
    ip <- (odataWSTr[1:(length(odataWSTr) - 1)] + odataWSTr[2:length(
        odataWSTr)])/2
    ip <- c(ip, (odataWSTr[length(odataWSTr)] + odataWSTr[1])/2)    #
#
# Now compute prediction error
#
    pe <- norm(ip, edata)   #
#
# Now repeat all the above the other way around.
#
#
# Build even thresholded estimate, then, threshold and rebuild
#
    edataWST <- wst(edata, filter.number = filter.number, family = family)
    edataWSTt <- threshold.wst(edataWST, levels = levels, policy = "manual",
        value = threshold, verbose = thverb)
    if(InverseType == "average")
        edataWSTr <- AvBasis.wst(edataWSTt) #
    else if(InverseType == "minent") {
        edataNV <- MaNoVe(edataWSTt)
        cat("EVEN Node Vector\n")
        cat("---------------\n")
        print(edataNV)
        edataWSTr <- InvBasis.wst(edataWSTt, nv = edataNV)
    }
    else stop(paste("Unknown InverseType: ", InverseType))
    ip <- (edataWSTr[1:(length(edataWSTr) - 1)] + edataWSTr[2:length(
        edataWSTr)])/2
    ip <- c(ip, (edataWSTr[length(edataWSTr)] + edataWSTr[1])/2)    #
#
# Now compute prediction error
#
    pe <- (pe + norm(ip, odata))/2
    if(verbose != 0) {
        cat("For threshold value\n")
        print(threshold)
        cat("The pe estimate is ", pe, "\n")
    }
    pe
}
"HaarConcat"<-
function()
{
    x1 <- HaarMA(n = 128, order = 1)
    x2 <- HaarMA(n = 128, order = 2)
    x3 <- HaarMA(n = 128, order = 3)
    x4 <- HaarMA(n = 128, order = 4)
    c(x1, x2, x3, x4)
}
"HaarMA"<-
function(n, sd = 1, order = 5)
{
#
#   Generate Haar MA realization
#
#   n - number of observations; sd=variance of increments; order=MA order
# 
    z <- rnorm(n = n + (2^order) - 1, mean = 0, sd = sd)
    J <- order
    x <- rep(0, n)
    for(i in (2^J):(2^(J - 1) + 1))
        x <- x + z[i:(n + i - 1)]
    for(i in (2^(J - 1)):1)
        x <- x - z[i:(n + i - 1)]
    x <- x * 2^( - J/2)
    return(x)
}
"InvBasis"<-
function(...)
UseMethod("InvBasis")
"InvBasis.wp"<-
function(wp, nvwp, pktlist, verbose = FALSE, ...)
{
    nlev <- nlevelsWT(wp)
    if(missing(pktlist)) {
        pktlist <- print.nvwp(nvwp, printing = FALSE)
        if(nlev != nlevelsWT(nvwp)) {
            stop("The node vector you supplied cannot have arisen from the wavelet packet object you supplied as they have different numbers of levels"
                )
        }
    }
    lpkts <- length(pktlist$level)
    ndata <- 2^nlev
    cfvc <- rep(0, ndata)
    ixvc <- cfvc
    counter <- 0
    for(i in 1:lpkts) {
        lev <- pktlist$level[i]
        pkt <- pktlist$pkt[i]
        coefs <- getpacket(wp, level = lev, index = pkt)
        pklength <- 2^lev
        pkleftix <- pkt * pklength + 1
        pkrightix <- pkleftix + pklength - 1
        cfvc[pkleftix:pkrightix] <- coefs
        ixvc[pkleftix:pkrightix] <- counter
        if(verbose == TRUE) {
            cat("Level: ", lev, "\n")
            cat("Packet: ", pkt, "\n")
            cat("coefs: ")
            print(coefs)
            cat("---\n")
            cat("Packet length: ", pklength, "\n")
            cat("Packet left ix: ", pkleftix, "\n")
            cat("Packet right ix: ", pkrightix, "\n")
            cat("ixvc: ")
            print(ixvc)
            cat("---\n")
            cat("cfvc: ")
            print(cfvc)
            cat("---\n")
        }
        counter <- counter + 1
    }
    if(verbose == TRUE) {
        cat("SWEEPER Stage\n")
    }
    sweeper <- rle(ixvc)$lengths
    mx <- min(sweeper)
    while(mx < ndata) {
        ix <- ((1:length(sweeper))[sweeper == mx])[1]
        csweeper <- cumsum(c(1, sweeper))[1:length(sweeper)]
        lix <- sweeper[ix]
        rix <- sweeper[ix + 1]
        if(lix != rix)
            stop(paste(
                "wavethresh error: lix and rix are not the same. lix is ",
                lix, " rix is ", rix))
        if(verbose == TRUE) {
            cat("Sweeper: ")
            print(sweeper)
            cat("Cumsum Sweeper: ")
            print(csweeper)
            cat("At sweeper index position ", ix, "\n")
            cat("Left ix is ", lix, "\n")
            cat("Right ix is ", rix, "\n")
            cat("Corresponds to ", csweeper[ix], csweeper[ix + 1], 
                "\n")
        }
        cfixl <- csweeper[ix]
        cfixr <- csweeper[ix + 1]
        pklength <- lix
        c.in <- cfvc[cfixl:(cfixl + pklength - 1)]
        d.in <- cfvc[cfixr:(cfixr + pklength - 1)]
        c.out <- conbar(c.in, d.in, wp$filter)
        cfvc[cfixl:(cfixr + pklength - 1)] <- c.out
        sweeper <- sweeper[ - ix]
        sweeper[ix] <- rix + lix
        mx <- min(sweeper)
    }
    cfvc
}
"InvBasis.wst"<-
function(wst, nv, ...)
{
#
#
# Perform an inverse on wst given specification in nv
#
# indexlist is a list of packet indices for access into appropriate levels of
# wst, nrsteps will be the number of reconstruction steps
#
    pnv <- print.nv(nv, printing = FALSE)
    indexlist <- rev(pnv$indexlist)
    rvector <- pnv$rvector
    nrsteps <- length(indexlist)    #
#
#   blevel is the bottom level in the decomposition
#
    blevel <- nlevelsWT(nv) - nrsteps #
#
#   Now extract the data and put it all in a vector
#
    rdata <- getpacket(wst, level = blevel, index = indexlist[1], type = 
        "C")
    ldata <- length(rdata)
    D <- getpacket(wst, level = blevel, index = indexlist[1])
    rdata <- c(rdata, D)
    ldata <- c(ldata, length(D))
    for(i in 2:nrsteps) {
        D <- getpacket(wst, level = (blevel + i - 1), index = indexlist[
            i])
        rdata <- c(rdata, D)
        ldata <- c(ldata, length(D))
    }
    error <- 0
    invswr <- .C("wavepackrecon",
        rdata = as.double(rdata),
        ldata = as.integer(ldata),
        nrsteps = as.integer(nrsteps),
        rvector = as.integer(rvector),
        H = as.double(wst$filter$H),
        LengthH = as.integer(length(wst$filter$H)),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(invswr$error != 0)
        stop(paste("Error code was ", invswr$error, 
            " from wavepackrecon"))
    return(invswr$rdata)
}
"IsEarly"<-
function(x)
UseMethod("IsEarly")
"IsEarly.default"<-
function(x)
{
    return(FALSE)
}
"IsEarly.wd"<-
function(x)
{
    if(is.null(x$date))
        return(TRUE)
    else return(FALSE)
}
"IsPowerOfTwo"<-
function(n)
{
    tvec <- (n == trunc(n))
    r <- log(n)/log(2)
    tvec <- tvec & (r == trunc(r))
    r[tvec == FALSE] <- NA
    r
}
"LocalSpec"<-
function(...)
UseMethod("LocalSpec")
"LocalSpec.wd"<-
function(wdS, lsmooth = "none", nlsmooth = FALSE, prefilter = TRUE, verbose = FALSE, 
    lw.number = wdS$filter$filter.number, lw.family = wdS$filter$family, 
    nlw.number = wdS$filter$filter.number, nlw.family = wdS$filter$family, 
    nlw.policy = "LSuniversal", nlw.levels = 0:(nlevelsWT(wdS) - 1), nlw.type
     = "hard", nlw.by.level = FALSE, nlw.value = 0, nlw.dev = var, nlw.boundary
     = FALSE, nlw.verbose = FALSE, nlw.cvtol = 0.01, nlw.Q = 0.050000000000000003, 
    nlw.alpha = 0.050000000000000003, nlw.transform = I, nlw.inverse = I, 
    debug.spectrum = FALSE, ...)
{
#
#
#   Check the class of the object
#
    cwd <- class(wdS)
    if(is.null(cwd) || cwd != "wd")
        stop("Object must be of class wd to perform energy computation"
            )
    else if(wdS$type != "station")
        stop("swd type should be station (nondecimated)")
    lnlevels <- nlevelsWT(wdS)
    N <- 2^lnlevels
    if(verbose == TRUE) cat("Original data length was:", N, "\n")   #
#
#   Decide whether to do no smoothing, Fourier smoothing or wavelet
#   linear smoothing.
#
    if(lsmooth == "none") {
#
#
#       Just square the coefficients in the wdS object
#
        if(verbose == TRUE) cat("Squaring coefficients on level: ")
        for(i in (lnlevels - 1):0) {
            if(verbose == TRUE)
                cat(i, " ")
            v <- accessD(wdS, level = i)
            if(debug.spectrum == TRUE)
                spectrum(v, spans = c(11, 9, 7))
            v <- v^2
            if(debug.spectrum == TRUE)
                spectrum(v, spans = c(11, 9, 7))
            wdS <- putD(wdS, level = i, v = v)
        }
        if(verbose == TRUE)
            cat("\n")
    }
    else if(lsmooth == "Fourier") {
#
#   Perform smoothing using Fourier methods.
#   For each level take the real cts Fourier transform and smooth
#   by removing a proportion of the coefficients and inverting the
#   transform.
#
#   The amount of smoothing is controlled by the fracsmooth variable
#   Initially this is set to 1/2 as the frequencies we want to remove
#   are 1/2 to 1. When we move a level up the frequencies we want to
#   remove are above 1/4 and so on. Note that smoothing starts at
#   level J-2 (not J-1 as these are the frequencies between 1 and 2
#   and I'm not sure what to do with these yet). 
#
#
        if(verbose == TRUE) {
            cat("Performing Fourier linear smoothing\n")
            cat("Processing level: ")
        }
        fracsmooth <- 1/2
        for(i in (lnlevels - 2):0) {
            v <- accessD(wdS, level = i)
            if(debug.spectrum == TRUE)
                spectrum(v, spans = c(11, 9, 7))
            if(verbose == TRUE) cat(i, " ") #
#
#   Do prefiltering if necessary. This low-passes the actual coefficients
#   to that the cut-off is at the highest frequency of the current
#   (Littlewood-Paley) wavelet.
#
            if(prefilter == TRUE) {
                if(verbose == TRUE)
                  cat("prefilter\n")
                vfft <- rfft(v)
                n <- length(vfft)
                start <- 1 + n * fracsmooth
                if(start <= n)
                  vfft[max(1, start):n] <- 0
                v <- rfftinv(vfft)
                if(debug.spectrum == TRUE)
                  spectrum(v, spans = c(11, 9, 7))
            }
#
#
#   Square the coefficients!
#
            v <- v^2
            if(debug.spectrum == TRUE) spectrum(v, spans = c(11, 9, 7)
                  ) #
#
#   Now carry out the Fourier smoothing. 
#
            vfft <- rfft(v)
            n <- length(vfft)
            start <- 1 + n * fracsmooth #
#
#               Maybe use something like this to adapt to
#               the shape of the wavelet?
#               start <- start * 0.77
#
            if(start <= n)
                vfft[max(1, start):n] <- 0
            v <- rfftinv(vfft)
            fracsmooth <- fracsmooth/2
            if(debug.spectrum == TRUE && i != 0)
                spectrum(v, spans = c(11, 9, 7))
            wdS <- putD(wdS, level = i, v = v)
        }
        if(verbose == TRUE)
            cat("\nSquaring top level only\n")
        v <- accessD(wdS, level = lnlevels - 1)
        if(debug.spectrum == TRUE)
            spectrum(v, spans = c(11, 9, 7))
        v <- v^2
        if(debug.spectrum == TRUE)
            spectrum(v, spans = c(11, 9, 7))
        wdS <- putD(wdS, level = lnlevels - 1, v)
    }
    else if(lsmooth == "wavelet") {
#
#
#   Do LINEAR wavelet smoothing
#
        if(verbose == TRUE) {
            cat("Performing LINEAR wavelet smoothing\n")
            cat("Processing level ")
        }
        fracsmooth <- 1/2
        for(i in 0:(lnlevels - 2)) {
            if(verbose == TRUE)
                cat(i, " ")
            v <- accessD(wdS, level = i)    #
#
#   Do prefiltering if necessary. This low-passes the actual coefficients
#   to that the cut-off is at the highest frequency of the current
#   (Littlewood-Paley) wavelet.
#
            if(debug.spectrum == TRUE)
                spectrum(v, spans = c(11, 9, 7))
            if(prefilter == TRUE) {
                if(verbose == TRUE)
                  cat("prefilter\n")
                vfft <- rfft(v)
                n <- length(vfft)
                start <- 1 + n * fracsmooth
                if(start <= n)
                  vfft[max(1, start):n] <- 0
                v <- rfftinv(vfft)
                if(debug.spectrum == TRUE)
                  spectrum(v, spans = c(11, 9, 7))
            }
#
#
#   Square the coefficients
#
            v <- v^2    #
#
#   Now do the linear wavelet smoothing. This takes each level (i), applies
#   the standard discrete wavelet transform and nulls levels higher than
#   the one we are at (j>i). The inverse transform is then applied and
#   the coefficients restored in the wdS object.
#
            if(debug.spectrum == TRUE)
                spectrum(v, spans = c(11, 9, 7))
            realwd <- wd(v, filter.number = lw.number, family = lw.family)
            realwd <- nullevels(realwd, levels = (i + 1):(nlevelsWT(
                realwd) - 1))
            v <- wr(realwd)
            if(debug.spectrum == TRUE && i != 0)
                spectrum(v, spans = c(11, 9, 7))
            wdS <- putD(wdS, level = i, v = v)
        }
        if(verbose == TRUE)
            cat("\nSquaring top level only\n")
        v <- accessD(wdS, level = lnlevels - 1)
        if(debug.spectrum == TRUE)
            spectrum(v, spans = c(11, 9, 7))
        v <- v^2
        if(debug.spectrum == TRUE)
            spectrum(v, spans = c(11, 9, 7))
        wdS <- putD(wdS, level = lnlevels - 1, v)
    }
    else stop(paste("Unknown lsmooth:", lsmooth))   #
    if(nlsmooth == TRUE) {
        if(verbose == TRUE) {
            cat("Performing non-linear wavelet smoothing\n")
            cat("Processing level: ")
        }
        for(i in ((lnlevels - 1):0)) {
            if(verbose == TRUE)
                cat(i, " ")
            v <- accessD(wdS, level = i)
            v <- nlw.transform(v)
            vwd <- wd(v, filter.number = nlw.number, family = nlw.family)
            vwdt <- threshold(vwd, levels = nlw.levels, type = 
                nlw.type, policy = nlw.policy, by.level = 
                nlw.by.level, value = nlw.value, dev = nlw.dev, 
                boundary = nlw.boundary, verbose = nlw.verbose, 
                cvtol = nlw.cvtol, Q = nlw.Q, alpha = nlw.alpha
                )
            v <- wr(vwdt)
            v <- nlw.inverse(v)
            wdS <- putD(wdS, level = i, v = v)
        }
        if(verbose == TRUE)
            cat("\n")
    }
    wdS
}
"LocalSpec.wst"<-
function(wst, ...)
{
    LocalSpec.wd(convert.wst(wst), ...)
}
"MaNoVe"<-
function(...)
UseMethod("MaNoVe")
"MaNoVe.wp"<-
function(wp, verbose = FALSE, ...)
{
    nlevels <- nlevelsWT(wp)
    LengthData <- dim(wp$wp)[[2]]
    upperctrl <- rep(0, LengthData - 1)
    upperl <- upperctrl
    firstl <- rev(c(0, cumsum(2^(0:(nlevels - 2)))))
    if(verbose == TRUE)
        verbose <- 1
    error <- 0
    tmp <- .C("wpCmnv",
        wp = as.double(wp$wp),
        LengthData = as.integer(LengthData),
        nlevels = as.integer(nlevels),
        upperctrl = as.integer(upperctrl),
        upperl = as.double(upperl),
        firstl = as.integer(firstl),
        verbose = as.integer(verbose),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(tmp$error != 0)
        stop(paste("Error condition ", tmp$error, 
            " reported from wpCmnv"))   #
    node.list <- vector("list", nlevels)
    matchcodes <- c("T", "B")
    vlength <- 2^(nlevels - 1)  #
#
#   Convert C to S
#
    firstl <- firstl + 1
    for(i in 1:nlevels) {
        first <- firstl[i]
        sv <- first:(first + vlength - 1)
        node.list[[i]]$upperctrl <- matchcodes[tmp$upperctrl[sv]]
        node.list[[i]]$upperl <- tmp$upperl[sv]
        vlength <- vlength/2
    }
    node.vector <- list(node.list = node.list, nlevels = nlevels)
    class(node.vector) <- "nvwp"
    node.vector
}
"MaNoVe.wst"<-
function(wst, entropy = Shannon.entropy, verbose = FALSE, stopper = FALSE, alg = "C", ...)
{
#
# Make a node vector. Use C code rather than the slow S code
#
    if(alg == "C") {
        if(verbose == TRUE)
            cat("Using C code version\n")
        nlevels <- nlevelsWT(wst) 
    #       node.vector <- vector("list", nlevels)
#       matchcodes <- c("S", "L", "R")
        LengthData <- dim(wst$wp)[[2]]
        upperctrl <- rep(0, LengthData - 1)
        upperl <- upperctrl
        firstl <- rev(c(0, cumsum(2^(0:(nlevels - 2)))))
        if(verbose == TRUE)
            verbose <- 1
        error <- 0
        tmp <- .C("Cmnv",
            wst = as.double(wst$wp),
            wstC = as.double(wst$Carray),
            LengthData = as.integer(LengthData),
            nlevels = as.integer(nlevels),
            upperctrl = as.integer(upperctrl),
            upperl = as.double(upperl),
            firstl = as.integer(firstl),
            verbose = as.integer(verbose),
            error = as.integer(error), PACKAGE = "wavethresh")
        if(tmp$error != 0)
            stop(paste("Error condition ", tmp$error, 
                " reported from Cmnv")) #
        node.list <- vector("list", nlevels)
        matchcodes <- c("S", "L", "R")
        vlength <- 2^(nlevels - 1)  #
#
#   Convert C to S
#
        firstl <- firstl + 1
        for(i in 1:nlevels) {
            first <- firstl[i]
            sv <- first:(first + vlength - 1)
            node.list[[i]]$upperctrl <- matchcodes[tmp$upperctrl[sv
                ]]
            node.list[[i]]$upperl <- tmp$upperl[sv]
            vlength <- vlength/2
        }
        node.vector <- list(node.list = node.list, nlevels = nlevels)
    }
    else {
        if(verbose == TRUE)
            cat("Using S code version\n")
        nlevels <- nlevelsWT(wst)
        node.vector <- vector("list", nlevels)
        matchcodes <- c("S", "L", "R")
        for(i in 0:(nlevels - 1)) {
            if(verbose == TRUE)
                cat("Lower level: ", i, "\n")
            nll <- 2^(nlevels - i)
            lowerl <- rep(0, nll)
            nul <- nll/2
            upperl <- rep(0, nul)
            upperctrl <- rep("", nul)
            if(verbose == TRUE)
                cat("Packets. Lower: ", nll, " Upper ", nul, 
                  "\n")
            for(j in 0:(nul - 1)) {
                if(verbose == TRUE)
                  cat("Upper level index: ", j, "\n")
                kl <- 2 * j
                kr <- 2 * j + 1
                mother.entropy <- entropy(getpacket(wst, level
                   = i + 1, index = j, type = "C"))
                if(i == 0) {
                  daughter.left.entropy <- entropy(c(getpacket(
                    wst, level = i, index = kl), getpacket(wst, 
                    level = i, index = kl, type = "C")))
                  daughter.right.entropy <- entropy(c(getpacket(
                    wst, level = i, index = kr), getpacket(wst, 
                    level = i, index = kr, type = "C")))
                }
                else {
                  if(verbose == TRUE)
                    cat("Left Ent C contrib ", node.vector[[i]]$
                      upperl[kl + 1], "\n")
                  daughter.left.entropy <- entropy(getpacket(
                    wst, level = i, index = kl)) + node.vector[[
                    i]]$upperl[kl + 1]
                  if(verbose == TRUE)
                    cat("Right Ent C contrib ", node.vector[[i
                      ]]$upperl[kr + 1], "\n")
                  daughter.right.entropy <- entropy(getpacket(
                    wst, level = i, index = kr)) + node.vector[[
                    i]]$upperl[kr + 1]
                }
                if(verbose == TRUE) {
                  cat("\tMother ent.:  ", mother.entropy, "\n")
                  cat("\tDaug. l .ent: ", daughter.left.entropy,
                    "\n")
                  cat("\tDaug. r .ent: ", 
                    daughter.right.entropy, "\n")
                }
                ents <- c(mother.entropy, daughter.left.entropy,
                  daughter.right.entropy)
                pos <- match(min(ents), ents)
                upperctrl[j + 1] <- matchcodes[pos]
                upperl[j + 1] <- min(ents)
                if(verbose == TRUE)
                  cat("\tSelected ", upperctrl[j + 1], upperl[j +
                    1], "\n")
                if(stopper == TRUE)
                  scan()
            }
            node.vector[[i + 1]] <- list(upperctrl = upperctrl, 
                upperl = upperl)
            if(verbose == TRUE)
                print(node.vector)
        }
        node.vector <- list(node.list = node.vector, nlevels = nlevels)
    }
    class(node.vector) <- "nv"
    node.vector
}
"PsiJ"<-
function(J, filter.number = 10, family = "DaubLeAsymm", tol = 1e-100, OPLENGTH
     = 10^7, verbose=FALSE)
{
    if (verbose==TRUE)
	    cat("Computing PsiJ\n")
    now <- proc.time()[1:2]
    if(J >= 0)
        stop("J must be negative integer")
    if(J - round(J) != 0)
        stop("J must be an integer")
    Psiorig <- Psiname(J = J, filter.number = filter.number, family = 
        family) #
#
#   See if matrix already exists. If so, return it
#
    if(exists(Psiorig, envir=WTEnv)) {
	if (verbose==TRUE)
		cat("Returning precomputed version\n")
        speed <- proc.time()[1:2] - now
	if (verbose==TRUE)
		cat("Took ", sum(speed), " seconds\n")
        return(get(Psiorig, envir=WTEnv))
    }
    H <- filter.select(filter.number = filter.number, family = family)$H
    wout <- rep(0, OPLENGTH)
    rlvec <- rep(0,  - J)
    error <- 0
    answer <- .C("PsiJ",
        J = as.integer( - J),
        H = as.double(H),
        LengthH = as.integer(length(H)),
        tol = as.double(tol),
        wout = as.double(wout),
        lwout = as.integer(length(wout)),
        rlvec = as.integer(rlvec),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(answer$error != 0) {
        if(answer$error == 160)
            cat("Increase ", OPLENGTH, " to be larger than ", 
                answer$lwout, "\n")
        stop(paste("Error code was ", answer$error))
    }
    speed <- proc.time()[1:2] - now
    if (verbose==TRUE)
	    cat("Took ", sum(speed), " seconds\n")
    m <- vector("list",  - J)
    lj <- c(0, cumsum(2 * answer$rlvec - 1))
    for(j in 1:( - J))
        m[[j]] <- answer$wout[(lj[j] + 1):lj[j + 1]]
    assign(Psiorig, m, envir=WTEnv)
    m
}
"PsiJmat"<-
function(J, filter.number = 10, family = "DaubLeAsymm", OPLENGTH = 10^7)
{
    J <-  - J
    P <- PsiJ( - J, filter.number = filter.number, family = family, 
        OPLENGTH = OPLENGTH)
    nc <- length(P[[J]])
    nr <- J
    m <- matrix(0, nrow = nr, ncol = nc)
    m[J,  ] <- P[[J]]
    for(j in 1:(J - 1)) {
        lj <- length(P[[j]])
        nz <- (nc - lj)/2
        z <- rep(0, nz)
        m[j,  ] <- c(z, P[[j]], z)
    }
    m
}
"Psiname"<-
function(J, filter.number, family)
{
    if(J >= 0)
        stop("J must be a negative integer")
    return(paste("Psi.",  - J, ".", filter.number, ".", family, sep = ""))
}
"ScalingFunction"<-
function(filter.number = 10, family = "DaubLeAsymm", resolution = 4096, 
    itlevels = 50)
{
    if(is.na(IsPowerOfTwo(resolution)))
        stop("Resolution must be a power of two")
    res <- 4 * resolution   #
#
# Select filter and work out some fixed constants
#
    H <- filter.select(filter.number = filter.number, family = family)$H
    lengthH <- length(H)
    ll <- lengthH
    v <- rep(0, res)    #
#
# Set initial coefficient to 1 in 2nd position on 1st level
#
    v[2] <- 1   #
#
# Now iterate the successive filtering operations to build up the scaling
# function. The actual filtering is carried out by the C routine CScalFn.
#
    for(it in 1:itlevels) {
        ans <- rep(0, res)
        z <- .C("CScalFn",
            v = as.double(v),
            ans = as.double(ans),
            res = as.integer(res),
            H = as.double(H),
            lengthH = as.integer(lengthH), PACKAGE = "wavethresh")  #
#
#       We only ever take the first half of the result
#
        v <- z$ans[1:(res/2)]   #
#
#       Set all other coefficients equal to zero. (This is because
#       rounding errors sometimes cause small values to appear).
#
        v[ - ((2^it + 1):(2^it + ll))] <- 0 
    #       plot(seq(from = 0, to = 2 * filter.number - 1, length = ll), v[(
#           2^it + 1):(2^it + ll)], type = "l")
        v <- sqrt(2) * v
        llbef <- ll
        vbef <- v   #
#
#       Check to see if the next iteration would send the number
#       of coefficients over the resolution that we can have.
#       Exit the loop if it does.
#
        if(2^(it + 1) + lengthH + ll * 2 - 2 > res/2) {
            cit <- it
            break
        }
#
#
#       ll is the number of coefficients that are nonzero in
#       any particular run. This formula updates ll for next time
#       round.
#
        ll <- lengthH + ll * 2 - 2  #
#
#       Add some zeroes to v to make it the right length.
#
        v <- c(v, rep(0, res - length(v)))
    }
    list(x = seq(from = 0, to = 2 * filter.number - 1, length = llbef), y
         = vbef[(2^cit + 1):(2^cit + llbef)])
}
"Shannon.entropy"<-
function(v, zilchtol = 1e-300)
{
    vsq <- v^2
    if(sum(vsq) < zilchtol)
        return(0)
    else {
        vsq[vsq == 0] <- 1
        return( - sum(vsq * log(vsq)))
    }
}
"TOgetthrda1"<-
function(dat, alpha)
{
    datsq <- sort(dat^2)
    a <- TOonebyone1(datsq, alpha)
    if(length(a) == length(datsq))
        if(1 - pchisq(datsq[1], 1) < alpha)
            ggg <- 0
        else ggg <- sqrt(datsq[1])
    else ggg <- sqrt(datsq[length(datsq) - length(a) + 1])
    return(ggg)
}
"TOgetthrda2"<-
function(dat, alpha)
{
    a <- TOonebyone2(dat, alpha)
    if(length(a) == length(dat))
        if(1 - pchisq(min(dat), 1) < alpha)
            ggg <- 0
        else ggg <- sqrt(min(dat))
    else ggg <- sqrt(max(dat[sort(order(dat)[1:(length(dat) - length(a) + 1
            )])]))
    return(ggg)
}
"TOkolsmi.chi2"<-
function(dat)
{
    n <- length(dat)
    return(max(abs(cumsum(dat) - ((1:n) * sum(dat))/n))/sqrt(2 * n))
}
"TOonebyone1"<-
function(dat, alpha)
{
    i <- length(dat)
    cc <- 1 - pchisq(dat[i], 1)^i
    while(cc[length(cc)] < alpha && i > 1) {
        i <- i - 1
        cc <- c(cc, 1 - pchisq(dat[i], 1)^i)
    }
    return(cc)
}
"TOonebyone2"<-
function(dat, alpha)
{
    crit <- c(seq(0.28000000000000003, 1.49, by = 0.01), seq(1.5, 2.48, by
         = 0.02))
    alph <- c(0.99999899999999997, 0.999996, 0.99999099999999996, 
        0.99997899999999995, 0.99995400000000001, 0.99990900000000005, 
        0.99982899999999997, 0.99969699999999995, 0.99948899999999996, 
        0.99917400000000001, 0.99871500000000002, 0.99807100000000004, 
        0.99719199999999997, 0.99602800000000002, 0.99452399999999996, 
        0.99262300000000003, 0.99026999999999998, 0.98741000000000001, 
        0.98399499999999995, 0.97997800000000002, 0.97531800000000002, 
        0.96998300000000004, 0.96394500000000005, 0.95718599999999998, 
        0.94969400000000004, 0.94146600000000003, 0.93250299999999997, 
        0.922817, 0.91242299999999998, 0.90134400000000003, 
        0.88960499999999998, 0.87724000000000002, 0.86428199999999999, 
        0.85077100000000005, 0.83677500000000005, 0.82224699999999995, 
        0.80732300000000001, 0.79201299999999997, 0.77636300000000003, 
        0.76041800000000004, 0.74421999999999999, 0.72781099999999999, 
        0.71123499999999995, 0.69452899999999995, 0.67773499999999998, 
        0.660887, 0.64401900000000001, 0.62716700000000003, 
        0.61036000000000001, 0.59362800000000004, 0.57699800000000001, 
        0.56049499999999997, 0.54414300000000004, 0.52795899999999996, 
        0.51197000000000004, 0.49619200000000002, 0.48063400000000001, 
        0.46531800000000001, 0.45025599999999999, 0.43545400000000001, 
        0.42093000000000003, 0.40668399999999999, 0.39273000000000002, 
        0.37907200000000002, 0.36571399999999998, 0.35266199999999998, 
        0.339918, 0.327484, 0.31536399999999998, 0.30355599999999999, 
        0.29205999999999999, 0.28087400000000001, 0.27000000000000002, 
        0.259434, 0.24917400000000001, 0.23921999999999999, 
        0.22956599999999999, 0.22020600000000001, 0.21113999999999999, 
        0.20236399999999999, 0.19387199999999999, 0.18565799999999999, 
        0.17771799999999999, 0.17005000000000001, 0.16264400000000001, 
        0.155498, 0.14860599999999999, 0.141962, 0.13555800000000001, 
        0.129388, 0.12345200000000001, 0.117742, 0.11225, 0.10697, 
        0.101896, 0.097028000000000003, 0.092352000000000004, 
        0.087868000000000002, 0.083568000000000003, 
        0.079444000000000001, 0.075495000000000007, 
        0.071711999999999998, 0.068092, 0.064630000000000007, 
        0.061317999999999998, 0.058152000000000002, 
        0.055128000000000003, 0.052243999999999999, 
        0.049487999999999997, 0.046857999999999997, 
        0.044350000000000001, 0.041959999999999997, 
        0.039682000000000002, 0.037513999999999999, 0.035448, 0.033484, 
        0.031618, 0.029842, 0.028153999999999998, 0.026551999999999999, 
        0.02503, 0.023588000000000001, 0.022218000000000002, 
        0.019689999999999999, 0.017422, 0.015389999999999999, 
        0.013573999999999999, 0.011952000000000001, 0.010508, 
        0.0092230000000000003, 0.0080829999999999999, 
        0.0070720000000000002, 0.0061770000000000002, 
        0.0053880000000000004, 0.0046909999999999999, 0.004078, 
        0.0035400000000000002, 0.003068, 0.0026540000000000001, 
        0.0022929999999999999, 0.001977, 0.0017030000000000001, 
        0.001464, 0.001256, 0.0010759999999999999, 
        0.00092100000000000005, 0.00078700000000000005, 
        0.00067100000000000005, 0.00057200000000000003, 0.000484, 
        0.00041199999999999999, 0.00035, 0.00029500000000000001, 
        0.00025000000000000001, 0.00021000000000000001, 
        0.00017799999999999999, 0.00014799999999999999, 0.000126, 
        0.00010399999999999999, 8.7999999999999998e-05, 
        7.3999999999999996e-05, 6.0000000000000002e-05, 5.1e-05, 
        4.1999999999999998e-05, 3.4999999999999997e-05, 
        3.0000000000000001e-05, 2.4000000000000001e-05, 
        2.0000000000000002e-05, 1.5999999999999999e-05, 
        1.2999999999999999e-05, 1.1e-05, 9.0000000000000002e-06)
    if(alpha < min(alph) || alpha > max(alph))
        stop("alpha =", alpha, "is out of range")
    ind <- match(TRUE, alpha > alph)
    critval <- crit[ind - 1] + ((alph[ind - 1] - alpha) * (crit[ind] - crit[
        ind - 1]))/(alph[ind - 1] - alph[ind])
    i <- length(dat)
    cc <- TOkolsmi.chi2(dat)
    while(cc[length(cc)] > critval && i > 1) {
        i <- i - 1
        cc <- c(cc, TOkolsmi.chi2(dat[sort(order(dat)[1:i])]))
    }
    return(cc)
}
"TOshrinkit"<-
function(coeffs, thresh)
{
    sign(coeffs) * pmax(abs(coeffs) - thresh, 0)
}
"TOthreshda1"<-
function(ywd, alpha = 0.050000000000000003, verbose = FALSE, return.threshold = FALSE)
{
    if(verbose)
        cat("Argument checking\n")
    ctmp <- class(ywd)
    if(is.null(ctmp))
        stop("ywd has no class")
    else if(ctmp != "wd")
        stop("ywd is not of class wd")
    if(alpha <= 0 || alpha >= 1)
        stop("alpha out of range")
    ans <- ywd
    n <- length(ywd$D)
    nlev <- log(n + 1, base = 2) - 1
    i <- nlev
    iloc <- 1
    while(i >= 0) {
        gg <- ywd$D[iloc:(iloc + 2^i - 1)]
        thresh <- TOgetthrda1(gg, alpha)
        if(verbose) {
            cat(paste("At level ", i, ", the threshold is ", thresh,
                "\n", sep = ""))
        }
        if(return.threshold)
            if(i == nlev)
                rt <- thresh
            else rt <- c(thresh, rt)
        else ans$D[iloc:(iloc + 2^i - 1)] <- TOshrinkit(ywd$D[iloc:(
                iloc + 2^i - 1)], thresh)
        iloc <- iloc + 2^i
        i <- i - 1
    }
    if(return.threshold)
        return(rt)
    else return(ans)
}
"TOthreshda2"<-
function(ywd, alpha = 0.050000000000000003, verbose = FALSE, return.threshold = FALSE)
{
    if(verbose)
        cat("Argument checking\n")
    ctmp <- class(ywd)
    if(is.null(ctmp))
        stop("ywd has no class")
    else if(ctmp != "wd")
        stop("ywd is not of class wd")
    if(alpha <= 9.0000000000000002e-06 || alpha >= 0.99999899999999997)
        stop("alpha out of range")
    ans <- ywd
    n <- length(ywd$D)
    nlev <- log(n + 1, base = 2) - 1
    i <- nlev
    iloc <- 1
    while(i >= 0) {
        gg <- ywd$D[iloc:(iloc + 2^i - 1)]
        thresh <- TOgetthrda2(gg^2, alpha)
        if(verbose) {
            cat(paste("At level ", i, ", the threshold is ", thresh,
                "\n", sep = ""))
        }
        if(return.threshold)
            if(i == nlev)
                rt <- thresh
            else rt <- c(thresh, rt)
        else ans$D[iloc:(iloc + 2^i - 1)] <- TOshrinkit(ywd$D[iloc:(
                iloc + 2^i - 1)], thresh)
        iloc <- iloc + 2^i
        i <- i - 1
    }
    if(return.threshold)
        return(rt)
    else return(ans)
}
"WaveletCV"<-
function(ynoise, x = 1:length(ynoise), filter.number = 10, family = 
    "DaubLeAsymm", thresh.type = "soft", tol = 0.01, verbose = 0, plot.it
     = TRUE, ll = 3)
{
#
#   Switch on verbosity for function calls if necessary
#
    if(verbose == 2) CallsVerbose <- TRUE else CallsVerbose <- FALSE
    if(verbose == 1)
        cat("WaveletCV: Wavelet model building\nThinking ")
    n <- length(ynoise)
    ywd <- wd(ynoise, filter.number = filter.number, family = family, 
        verbose = CallsVerbose)
    univ.threshold <- threshold(ywd, type = thresh.type, return.threshold
         = TRUE, lev = ll:(nlevelsWT(ywd) - 1), verbose = CallsVerbose,
	policy="universal")[1]
    if(verbose == 1) {
        cat("Universal threshold: ", univ.threshold, "\n")
        cat("Now doing universal threshold reconstruction...")
    }
    yuvtwd <- threshold(ywd, type = thresh.type, lev = ll:(nlevelsWT(ywd) - 1),
        verbose = CallsVerbose, policy="universal")
    if(verbose == 1)
        cat("done\nNow reconstructing...")
    yuvtwr <- wr(yuvtwd, verbose = CallsVerbose)
    if(verbose == 1)
        cat("done\nNow plotting universal thresholded\n")
    if(plot.it == TRUE) {
        oldpar <- par(mfrow = c(2, 2))
        matplot(x, cbind(ynoise, yuvtwr), type = "l", main = 
            "Universal Threshold Reconstruction", xlab = "x", col
             = c(3, 2), lty = c(3, 2))
    }
    if(verbose == 1)
        cat("Now optimising cross-validated error estimate\n")
    R <- 0.61803399000000003
    C <- 1 - R
    ax <- 0
    bx <- univ.threshold/2
    cx <- univ.threshold
    x0 <- ax
    x3 <- cx
    if(abs(cx - bx) > abs(bx - ax)) {
        x1 <- bx
        x2 <- bx + C * (cx - bx)
    }
    else {
        x2 <- bx
        x1 <- bx - C * (bx - ax)
    }
    fa <- rsswav(ynoise, value = ax, filter.number = filter.number, family
         = family, thresh.type = thresh.type, ll = ll)$ssq
    fb <- rsswav(ynoise, value = bx, filter.number = filter.number, family
         = family, thresh.type = thresh.type, ll = ll)$ssq
    fc <- rsswav(ynoise, value = cx, filter.number = filter.number, family
         = family, thresh.type = thresh.type, ll = ll)$ssq
    f1 <- rsswav(ynoise, value = x1, filter.number = filter.number, family
         = family, thresh.type = thresh.type, ll = ll)$ssq
    f2 <- rsswav(ynoise, value = x2, filter.number = filter.number, family
         = family, thresh.type = thresh.type, ll = ll)$ssq
    xkeep <- c(ax, cx, x1, x2)
    fkeep <- c(fa, fc, f1, f2)
    if(plot.it == TRUE) {
        plot(c(ax, bx, cx), c(fa, fb, fc))
        text(c(x1, x2), c(f1, f2), lab = c("1", "2"))
    }
    cnt <- 3
    while(abs(x3 - x0) > tol * (abs(x1) + abs(x2))) {
        cat("x0=", x0, "x1=", x1, "x2=", x2, "x3=", x3, "\n")
        cat("f1=", f1, "f2=", f2, "\n")
        if(f2 < f1) {
            x0 <- x1
            x1 <- x2
            x2 <- R * x1 + C * x3
            f1 <- f2
            f2 <- rsswav(ynoise, value = x2, filter.number = 
                filter.number, family = family, thresh.type = 
                thresh.type, ll = ll)
            if(verbose == 2) {
                cat("SSQ: ", signif(f2$ssq, 3), " DF: ", f2$df, 
                  "\n")
            }
            else if(verbose == 1)
                cat(".")
            f2 <- f2$ssq
            xkeep <- c(xkeep, x2)
            fkeep <- c(fkeep, f2)
            if(plot.it == TRUE)
                text(x2, f2, lab = as.character(cnt))
            cnt <- cnt + 1
        }
        else {
            x3 <- x2
            x2 <- x1
            x1 <- R * x2 + C * x0
            f2 <- f1
            f1 <- rsswav(ynoise, value = x1, filter.number = 
                filter.number, family = family, thresh.type = 
                thresh.type, ll = ll)
            if(verbose == 2)
                cat("SSQ: ", signif(f1$ssq, 3), " DF: ", f1$df, 
                  "\n")
            else if(verbose == 1)
                cat(".")
            f1 <- f1$ssq
            xkeep <- c(xkeep, x1)
            fkeep <- c(fkeep, f1)
            if(plot.it == TRUE)
                text(x1, f1, lab = as.character(cnt))
            cnt <- cnt + 1
        }
    }
    if(f1 < f2)
        tmp <- x1
    else tmp <- x2
    x1 <- tmp/sqrt(1 - log(2)/log(n))
    if(verbose == 1)
        cat("Correcting to ", x1, "\n")
    else if(verbose == 1)
        cat("\n")
    xvwd <- threshold(ywd, policy = "manual", value = x1, type = 
        thresh.type, lev = ll:(nlevelsWT(ywd)- 1))
    xvwddof <- dof(xvwd)
    xvwr <- wr(xvwd)
    if(plot.it == TRUE)
        matplot(x, cbind(ynoise, yuvtwr, xvwr), type = "l", main = 
            "XV Threshold Reconstruction", xlab = "x", col = c(3, 2,
            1))
    g <- sort.list(xkeep)
    xkeep <- xkeep[g]
    fkeep <- fkeep[g]
    list(x = x, ynoise = ynoise, xvwr = xvwr, yuvtwr = yuvtwr, xvthresh = 
        x1, uvthresh = univ.threshold, xvdof = xvwddof, uvdof = dof(
        yuvtwd), xkeep = xkeep, fkeep = fkeep)
}
"Whistory"<-
function(...)
UseMethod("Whistory")
"Whistory.wst"<-
function(wst, all = FALSE, ...)
{
    ntimes <- length(wst$date)
    if(ntimes == 1)
        cat("This object has not been modified\n")
    cat("This object has been modified ", ntimes - 1, " times\n")
    cat("The date of the last mod was ", wst$date[ntimes], "\n")
    cat("That modification was\n")
    cat(wst$history[ntimes - 1], "\n")
    if(all == TRUE) {
        cat("Complete history\n")
        cat("Modification dates\n")
        for(i in 1:ntimes)
            cat(wst$date[i], "\n")
        cat("Modification record\n")
        for(i in 1:ntimes)
            cat(wst$history[i - 1], "\n")
    }
    invisible()
}
"accessC"<-
function(...)
UseMethod("accessC")
"accessC.mwd"<-
function(mwd, level = nlevelsWT(mwd), ...)
{
#
#  Get smoothed data from multiple wavelet structure.
#
    ctmp <- class(mwd)
    if(is.null(ctmp))
        stop("mwd has no class")
    else if(ctmp != "mwd")
        stop("mwd is not of class mwd")
    if(level < 0)
        stop("Must have a positive level")
    else if(level > nlevelsWT(mwd))
        stop("Cannot exceed maximum number of levels")
    level <- level + 1
    first.last.c <- mwd$fl.dbase$first.last.c
    first.level <- first.last.c[level, 1]
    last.level <- first.last.c[level, 2]
    offset.level <- first.last.c[level, 3]
    n <- last.level + 1 - first.level
    coeffs <- mwd$C[, (offset.level + 1):(offset.level + n)]
    return(coeffs)
}
"accessC.wd"<-
function(wd, level = nlevelsWT(wd), boundary = FALSE, aspect = "Identity", ...)
{
    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
    ctmp <- class(wd)
    if(is.null(ctmp))
        stop("wd has no class")
    else if(ctmp != "wd")
        stop("wd is not of class wd")
    if(level < 0)
        stop("Must have a positive level")
    else if(level > nlevelsWT(wd))
        stop(paste("Cannot exceed maximum number of levels", nlevelsWT(wd)
            ))
    if(wd$bc == "interval") {
        if(level != wd$current.scale)
            stop(paste(
                "Requested wd object was decomposed to level ", 
                wd$current.scale, 
                " and so for \"wavelets on the interval\" objects I can only show this level for the scaling function coefficients\n"
                ))
        first.level <- wd$fl.dbase$first.last.c[1]
        last.level <- wd$fl.dbase$first.last.c[2]
        offset.level <- wd$fl.dbase$first.last.c[3]
        n <- last.level - first.level + 1
        coefs <- wd$transformed.vector[(offset.level + 1 - first.level):
            (offset.level + n - first.level)]
    }
    else {
        level <- level + 1
        first.last.c <- wd$fl.dbase$first.last.c
        first.level <- first.last.c[level, 1]
        last.level <- first.last.c[level, 2]
        offset.level <- first.last.c[level, 3]
        if(boundary == TRUE) {
            n <- last.level - first.level + 1
            coefs <- wd$C[(offset.level + 1):(offset.level + n)]
        }
        else {
            type <- wd$type
            if(type == "wavelet")
                n <- 2^(level - 1)
            else if(type == "station")
                n <- 2^nlevelsWT(wd)
            else stop("Unknown type component")
            coefs <- wd$C[(offset.level + 1 - first.level):(
                offset.level + n - first.level)]
        }
    }
    if(!is.null(tsp(wd)))
        tsp(coefs) <- tsp(wd)
    if(aspect == "Identity")
        return(coefs)
    else {
        fn <- get(aspect)
        return(fn(coefs))
    }
}
"accessC.wp"<-
function(wp, ...)
{
    stop("A wavelet packet object does not have ``levels'' of father wavelet coefficients. Use accessD to obtain levels of father and mother coefficients"
        )
}
"accessC.wst"<-
function(wst, level, aspect = "Identity", ...)
{
#
#
# Get all coefficients at a particular level
# First work out how many packets there are at this level
#
    nlevels <- nlevelsWT(wst)
    if(level < 0)
        stop("level must nonnegative")
    else if(level > nlevels)
        stop(paste("level must be smaller than ", nlevels - 1))
    coefs <- wst$Carray[level + 1,  ]
    if(aspect == "Identity")
        return(coefs)
    else {
        fn <- get(aspect)
        return(fn(coefs))
    }
}
"accessD"<-
function(...)
UseMethod("accessD")
"accessD.mwd"<-
function(mwd, level, ...)
{
#
# Get wavelet coefficients from multiple wavelet structure
#
    ctmp <- class(mwd)
    if(is.null(ctmp))
        stop("mwd has no class")
    else if(ctmp != "mwd")
        stop("mwd is not of class mwd")
    if(level < 0)
        stop("Must have a positive level")
    else if(level > (nlevelsWT(mwd) - 1))
        stop("Cannot exceed maximum number of levels")
    level <- level + 1
    first.last.d <- mwd$fl.dbase$first.last.d
    first.level <- first.last.d[level, 1]
    last.level <- first.last.d[level, 2]
    offset.level <- first.last.d[level, 3]
    n <- last.level + 1 - first.level
    coeffs <- mwd$D[, (offset.level + 1):(offset.level + n)]
    return(coeffs)
}
"accessD.wd"<-
function(wd, level, boundary = FALSE, aspect = "Identity", ...)
{
    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
    ctmp <- class(wd)
    if(is.null(ctmp))
        stop("wd has no class")
    else if(ctmp != "wd")
        stop("wd is not of class wd")
    if(level < 0)
        stop("Must have a positive level")
    else if(level > (nlevelsWT(wd) - 1))
        stop(paste("Cannot exceed maximum number of levels: ", wd$
            nlevels - 1))
    if(wd$bc == "interval") {
        level <- level - wd$current.scale
        objname <- deparse(substitute(wd))
        if(level < 0)
            stop(paste("The wd object: ", objname, 
                " was only decomposed down to level: ", wd$
                current.scale, " Try a larger level"))
        if(boundary == TRUE)
            stop("There are no boundary elements in a wavelets on the interval transform!"
                )
    }
    level <- level + 1
    first.last.d <- wd$fl.dbase$first.last.d
    first.level <- first.last.d[level, 1]
    last.level <- first.last.d[level, 2]
    offset.level <- first.last.d[level, 3]
    if(boundary == TRUE) {
        n <- last.level - first.level + 1
        coefs <- wd$D[(offset.level + 1):(offset.level + n)]
    }
    else {
        type <- wd$type
        if(type == "wavelet") {
            n <- 2^(level - 1)
            if(wd$bc == "interval")
                n <- last.level - first.level + 1
        }
        else if(type == "station")
            n <- 2^nlevelsWT(wd)
        else stop("Unknown type component")
        if(wd$bc != "interval")
            coefs <- wd$D[(offset.level + 1 - first.level):(
                offset.level + n - first.level)]
        else coefs <- wd$transformed.vector[(offset.level + 1 - 
                first.level):(offset.level + n - first.level)]
    }
    if(!is.null(tsp(wd)))
        tsp(coefs) <- tsp(wd)
    if(aspect == "Identity")
        return(coefs)
    else {
        fn <- get(aspect)
        return(fn(coefs))
    }
}
"accessD.wd3D"<-
function(obj, level = nlevelsWT(obj) - 1, block, ...)
{
    if(level < 0)
        stop(paste("Level cannot be accessed. You tried to access level",
            level, ". The minimum is zero"))
    else if(level >= nlevelsWT(obj))
        stop(paste("Level cannot be accessed. You tried to access level",
            level, ". The maximum level is", nlevelsWT(obj) - 1))
    halfsize <- 2^level
    size <- dim(obj$a)[1]
    GHH <- HGH <- GGH <- HHG <- GHG <- HGG <- GGG <- array(0, dim = rep(
        halfsize, 3))
    answer <- .C("getARRel",
        Carray = as.double(obj$a),
        size = as.integer(size),
        level = as.integer(level),
        GHH = as.double(GHH),
        HGH = as.double(HGH),
        GGH = as.double(GGH),
        HHG = as.double(HHG),
        GHG = as.double(GHG),
        HGG = as.double(HGG),
        GGG = as.double(GGG), PACKAGE = "wavethresh")
    thedim <- rep(halfsize, 3)  #
#
# Return HHH if level = 0
#
    if(missing(block)) {
        if(level == 0)
            list(HHH = array(obj$a[1, 1, 1], dim = thedim), GHH = 
                array(answer$GHH, dim = thedim), HGH = array(
                answer$HGH, dim = thedim), GGH = array(answer$
                GGH, dim = thedim), HHG = array(answer$HHG, dim
                 = thedim), GHG = array(answer$GHG, dim = 
                thedim), HGG = array(answer$HGG, dim = thedim), 
                GGG = array(answer$GGG, dim = thedim))
        else list(GHH = array(answer$GHH, dim = thedim), HGH = array(
                answer$HGH, dim = thedim), GGH = array(answer$
                GGH, dim = thedim), HHG = array(answer$HHG, dim
                 = thedim), GHG = array(answer$GHG, dim = 
                thedim), HGG = array(answer$HGG, dim = thedim), 
                GGG = array(answer$GGG, dim = thedim))
    }
    else {
        if(level != 0 && block == "HHH")
            stop("HHH only exists at level 0")
        else return(switch(block,
                HHH = array(obj$a[1, 1, 1], dim = thedim),
                GHH = array(answer$GHH, dim = thedim),
                HGH = array(answer$HGH, dim = thedim),
                GGH = array(answer$GGH, dim = thedim),
                HHG = array(answer$HHG, dim = thedim),
                GHG = array(answer$GHG, dim = thedim),
                HGG = array(answer$HGG, dim = thedim),
                GGG = array(answer$GGG, dim = thedim)))
    }
}
"accessD.wp"<-
function(wp, level, ...)
{
#
#
# Get all coefficients at a particular level
# First work out how many packets there are at this level
#
    nlev <- nlevelsWT(wp)
    if(level < 0)
        stop("level must nonnegative")
    else if(level > nlev - 1)
        stop(paste("level must be smaller than ", nlev - 1))
    npx <- 2^(nlev - level)
    return(wp$wp[level + 1,  ])
}
"accessD.wpst"<-
function(wpst, level, index, ...)
{
    nlev <- nlevelsWT(wpst)
    if(level < 0)
        stop("Level must be greater than or equal to 0")
    else if(level >= nlev)
        stop(paste("Level must be less than ", nlev))
    nwppkt <- 2^(nlev - level)  #
#
#   Check that packet index "index" is in range
#
    if(index < 0)
        stop("index must be greater than or equal to 0")
    else if(index >= nwppkt)
        stop(paste("index must be less than ", nwppkt))
    primary.index <- c2to4(index)   #
#
#   Now compute extra multiples for lower levels
#
    for(i in level:(nlev - 1)) {
        em <- 2^(2 * nlev - 2 * i - 1)
        primary.index <- c(primary.index, em + primary.index)
    }
#
#
#   Prepare some room for the answer
#
    weave <- rep(0, 2^nlev)
    ans <- .C("accessDwpst",
        coefvec = as.double(wpst$wpst),
        lansvec = as.integer(length(wpst$wpst)),
        nlev = as.integer(nlev),
        avixstart = as.integer(wpst$avixstart),
        primary.index = as.integer(primary.index),
        nwppkt = as.integer(nwppkt),
        pklength = as.integer(2^level),
        level = as.integer(level),
        weave = as.double(weave),
        lweave = as.double(length(weave)),
        error = as.integer(0), PACKAGE = "wavethresh")
    ans$weave
}
"accessD.wst"<-
function(wst, level, aspect = "Identity", ...)
{
#
#
# Get all coefficients at a particular level
# First work out how many packets there are at this level
#
    nlevels <- nlevelsWT(wst)
    if(level < 0)
        stop("level must nonnegative")
    else if(level > nlevels - 1)
        stop(paste("level must be smaller than ", nlevels - 1))
    npx <- 2^(nlevels - level)
    coefs <- wst$wp[level + 1,  ]
    if(aspect == "Identity")
        return(coefs)
    else {
        fn <- get(aspect)
        return(fn(coefs))
    }
}
"accessc"<-
function(irregwd.structure, level, boundary = FALSE)
{
    ctmp <- class(irregwd.structure)
    if(is.null(ctmp))
        stop("irregwd.structure has no class")
    else if(ctmp != "irregwd")
        stop("irregwd.structure is not of class irregwd")
    if(level < 0)
        stop("Must have a positive level")
    else if(level > (nlevelsWT(irregwd.structure) - 1))
        stop("Cannot exceed maximum number of levels")
    level <- level + 1
    first.last.d <- irregwd.structure$fl.dbase$first.last.d
    first.level <- first.last.d[level, 1]
    last.level <- first.last.d[level, 2]
    offset.level <- first.last.d[level, 3]
    if(boundary == TRUE) {
        n <- last.level - first.level + 1
        coefs <- irregwd.structure$c[(offset.level + 1):(offset.level + 
            n)]
    }
    else {
        n <- 2^(level - 1)
        coefs <- irregwd.structure$c[(offset.level + 1 - first.level):(
            offset.level + n - first.level)]
    }
    return(coefs)
}
"addpkt"<-
function(level, index, density, col, yvals)
{
    if(density < 0 || density > 1)
        stop("Density should be between 0 and 1")
    density <- density * 40
    y <- level
    level <- level - 1
    pktlength <- 2^level
    x <- index * pktlength
    h <- 1
    w <- pktlength
    if(missing(yvals))
        drawbox(x, y, w, h, density = density, col = col)
    else {
        xco <- seq(from = x, to = x + w, length = length(yvals))
        yco <- y + h/2 + (h * yvals)/(2 * max(abs(yvals)))
        lines(xco, yco)
    }
}
"av.basis"<-
function(wst, level, ix1, ix2, filter)
{
    if(level != 0) {
        cl <- conbar(av.basis(wst, level - 1, 2 * ix1, 2 * ix1 + 1, 
            filter), getpacket(wst, level = level, index = ix1), 
            filter = filter)
        cr <- rotateback(conbar(av.basis(wst, level - 1, 2 * ix2, 2 * 
            ix2 + 1, filter), getpacket(wst, level = level, index
             = ix2), filter = filter))
    }
    else {
        cl <- conbar(getpacket(wst, level = level, index = ix1, type = 
            "C"), getpacket(wst, level = level, index = ix1), 
            filter)
        cr <- rotateback(conbar(getpacket(wst, level = level, index = 
            ix2, type = "C"), getpacket(wst, level = level, index
             = ix2), filter))
    }
    return(0.5 * (cl + cr))
}
"basisplot"<-
function(x, ...)
UseMethod("basisplot")
"basisplot.BP"<-
function(x, num = min(10, length(BP$level)), ...)
{
	BP <- x
    plotpkt(nlevelsWT(BP))
    dnsvec <- BP$basiscoef[1:num]
    dnsvec <- dnsvec/max(abs(dnsvec))
    for(i in 1:num)
        addpkt(BP$level[i], BP$pkt[i], dnsvec[i], col = 1)
}
"basisplot.wp"<-
function(x, draw.mode = FALSE, ...)
{
	wp <- x
    J <- nlevelsWT(wp)
    oldl <- -1
    zero <- rep(0, 2^J)
    rh <- 2^(J - 1)
    zwp <- wp(zero, filter.number = wp$filter$filter.number, family = wp$
        filter$family)
    plotpkt(J)
    for(j in 0:(J - 1))
        for(k in 0:(2^(J - j) - 1))
            addpkt(j, k, 0, col = 1)
    znv <- MaNoVe(zwp)
    origznv <- znv
    cat("Select packets: Left: select. Right: exit\n")
    endit <- 0
    while(endit == 0) {
        n <- locator(n = 1)
        if(length(n) == 0)
            endit <- 1
        else {
            sellevel <- floor(n$y)
            if(sellevel < 1 || sellevel > (J - 1))
                cat("Click on shaded boxes\n")
            else {
                npkts <- 2^(J - sellevel)
                if(n$x < 0 || n$x > rh)
                  cat("Click on shaded boxes\n")
                else {
                  pknumber <- floor((npkts * n$x)/rh)
                  if(draw.mode == TRUE && oldl > -1) {
                    addpkt(oldl, oldpn, 1, col = 3)
                  }
                  addpkt(sellevel, pknumber, 1, col = 2)
                  znv$node.list[[sellevel]]$upperctrl[pknumber + 
                    1] <- "T"
                  if(draw.mode == TRUE) {
                    oldl <- sellevel
                    oldpn <- pknumber
                    pktl <- 2^sellevel
                    nhalf <- floor(pktl/2)
                    pkt <- c(rep(0, nhalf), 1, rep(0, nhalf - 1
                      ))
                    nzwp <- putpacket(zwp, level = sellevel, 
                      index = pknumber, packet = pkt)
                    cat("Computing WAIT...")
                    ans <- InvBasis(nzwp, nv = znv)
                    cat("d o n e.\n")
                    znv <- origznv
                    dev.set()
                    ts.plot(ans, xlab = "x", ylab = 
                      "Wavelet packet basis function")
                    dev.set()
                  }
                }
            }
        }
    }
    znv
}
"c2to4"<-
function(index)
{
#
# Represent index in base 2. Then use this representation and think of
# it in base 4 to get the number
#
    ans <- .C("c2to4",
        index = as.integer(index),
        answer = as.integer(0) ,PACKAGE = "wavethresh")
    ans$answer
}
"compare.filters"<-
function(f1, f2)
{
    if(f1$family != f2$family)
        return(FALSE)
    else if(f1$filter.number != f2$filter.number)
        return(FALSE)
    else return(TRUE)
}
"compress"<-
function(...)
UseMethod("compress")
"compress.default"<-
function(v, verbose = FALSE, ...)
{
    n <- length(v)
    r <- sum(v != 0)
    if(n > 2 * r) {
        position <- (1:n)[v != 0]
        values <- v[position]
        answer <- list(position = position, values = values, 
            original.length = n)
        class(answer) <- "compressed"
        if(verbose == TRUE)
            cat("Compressed ", n, " into ", 2 * r, "(", signif((100 *
                2 * r)/n, 3), "%)\n")
        return(answer)
    }
    else {
        answer <- list(vector = v)
        class(answer) <- "uncompressed"
        if(verbose == TRUE)
            cat("No compression\n")
        return(answer)
    }
}
"compress.imwd"<-
function(x, verbose = FALSE, ...)
{
    if(verbose == TRUE) cat("Argument checking...") #
#
#       Check class of imwd
#
    if(verbose == TRUE)
        cat("Argument checking\n")
    ctmp <- class(x)
    if(is.null(ctmp))
        stop("imwd has no class")
    else if(ctmp != "imwd")
        stop("imwd is not of class imwd")
    squished <- list(nlevels = nlevelsWT(x), fl.dbase = x$fl.dbase, 
        filter = x$filter, w0Lconstant = x$w0Lconstant, type = 
        x$type, bc = x$bc)    #
#
#   Go round loop compressing each set of coefficients
#
    for(level in 0:(nlevelsWT(x) - 1)) {
        if(verbose == TRUE)
            cat("Level ", level, "\n\t")
        nm <- lt.to.name(level, "CD")
        if(verbose == TRUE)
            cat("CD\t")
        squished[[nm]] <- compress.default(x[[nm]], verbose = verbose)
        nm <- lt.to.name(level, "DC")
        if(verbose == TRUE)
            cat("\tDC\t")
        squished[[nm]] <- compress.default(x[[nm]], verbose = verbose)
        nm <- lt.to.name(level, "DD")
        if(verbose == TRUE)
            cat("\tDD\t")
        squished[[nm]] <- compress.default(x[[nm]], verbose = verbose)
    }
    class(squished) <- c("imwdc")
    if(verbose == TRUE)
        cat("Overall compression: Was: ", w <- object.size(x), 
            " Now:", s <- object.size(squished), " (", signif((100 * 
            s)/w, 3), "%)\n")
    squished
}
"conbar"<-
function(c.in, d.in, filter)
{
#
# S interface to C routine conbar
#
    LengthCout <- 2 * length(c.in)
    c.out <- rep(0, LengthCout)
    answer <- .C("conbarL",
        c.in = as.double(c.in),
        LengthCin = as.integer(length(c.in)),
        firstCin = as.integer(0),
        d.in = as.double(d.in),
        LengthDin = as.integer(length(d.in)),
        firstDin = as.integer(0),
        H = as.double(filter$H),
        LengthH = as.integer(length(filter$H)),
        c.out = as.double(c.out),
        LengthCout = as.integer(LengthCout),
        firstCout = as.integer(0),
        lastCout = as.integer(LengthCout - 1),
        type = as.integer(1),
        bc = as.integer(1), PACKAGE = "wavethresh")
    answer$c.out
}
"convert"<-
function(...)
UseMethod("convert")
"convert.wd"<-
function(wd, ...)
{
#
#
# Convert a wd station object into a wst object
#
#
# First create object of same size and type of desired return object.
#
    if(wd$type != "station") stop(
            "Object to convert must be of type \"station\" ")
    n <- 2^nlevelsWT(wd)
    dummy <- rep(0, n)
    tmpwst <- wst(dummy, filter.number = wd$filter$filter.number, family = wd$
        filter$family)
    tmpwst$date <- wd$date  #
#
#   Now we've got the skeleton let's fill in all the details.
#
    arrvec <- getarrvec(nlevelsWT(wd), sort = FALSE)
    for(lev in (nlevelsWT(wd) - 1):1) {
        ds <- accessD.wd(wd, level = lev)
        cs <- accessC.wd(wd, level = lev)
        ds <- ds[arrvec[, nlevelsWT(wd) - lev]]
        cs <- cs[arrvec[, nlevelsWT(wd) - lev]]
        tmpwst <- putD(tmpwst, level = lev, v = ds)
        tmpwst <- putC(tmpwst, level = lev, v = cs)
    }
#
#
#   And put final level in for Cs and Ds (for wst only)
#
    tmpwst <- putC(tmpwst, level = nlevelsWT(wd), v = accessC(wd, level = wd$
        nlevels))   #
    tmpwst <- putD(tmpwst, level = nlevelsWT(wd), v = accessC(wd, level = wd$
        nlevels))   #
#
#   And zeroth level
#
    tmpwst <- putC(tmpwst, level = 0, v = accessC(wd, level = 0))
    arrvec <- sort.list(levarr(1:n, levstodo = nlevelsWT(wd)))
    tmpwst <- putD(tmpwst, level = 0, v = accessD(wd, level = 0)[arrvec])
    tmpwst
}
"convert.wst"<-
function(wst, ...)
{
#
#
# Convert a wst object into a wd type station object
#
#
# First create object of same size and type of desired return object.
#
    n <- 2^nlevelsWT(wst)
    dummy <- rep(0, n)
    tmpwd <- wd(dummy, type = "station", filter.number = wst$filter$filter.number, 
        family = wst$filter$family)
    tmpwd$date <- wst$date  #
#
#   Now we've got the skeleton let's fill in all the details.
#
    arrvec <- getarrvec(nlevelsWT(wst))
    for(lev in (nlevelsWT(wst) - 1):1) {
        ds <- accessD.wst(wst, level = lev)
        cs <- accessC.wst(wst, level = lev)
        ds <- ds[arrvec[, nlevelsWT(wst) - lev]]
        cs <- cs[arrvec[, nlevelsWT(wst) - lev]]
        ixs <- putD(tmpwd, level = lev, v = ds, index = TRUE)
        tmpwd$D[ixs$ix1:ixs$ix2] <- ds
        ixs <- putC(tmpwd, level = lev, v = cs, index = TRUE)
        tmpwd$C[ixs$ix1:ixs$ix2] <- cs
    }
#
#
#   And put final level in for Cs
#
    tmpwd <- putC(tmpwd, level = nlevelsWT(wst), v = accessC(wst, level = wst$
        nlevels))   #
#
#   And zeroth level
#
    tmpwd <- putC(tmpwd, level = 0, v = accessC(wst, level = 0))
    arrvec <- levarr(1:n, levstodo = nlevelsWT(wst))
    tmpwd <- putD(tmpwd, level = 0, v = accessD(wst, level = 0)[arrvec])
    tmpwd
}
"dof"<-
function(wd)
{
    cwd <- class(wd)
    if(is.null(cwd)) {
        stop("Object has no class")
    }
    else if(cwd != "wd")
        stop("Object is not of class wd")
    else {
#
# Count number of non-zero coefficients
#
        nlev <- nlevelsWT(wd) #
#
#   nnonzero counts the number of nonzero coefficients
#   This is already 1, since the C contains first level constant
#
        nnonzero <- 1
        for(i in 0:(nlev - 1)) {
            nnonzero <- nnonzero + sum(accessD(wd, lev = i) != 0)
        }
    }
    nnonzero
}
"doppler"<-
function(t)
{
    sqrt(t * (1 - t)) * sin((2 * pi * 1.05)/(t + 0.050000000000000003))
}
"draw"<-
function(...)
UseMethod("draw")
"draw.default"<-
function(filter.number = 10, family = "DaubLeAsymm", resolution = 8192, verbose
     = FALSE, plot.it = TRUE, main = "Wavelet Picture", sub = zwd$filter$name, 
    xlab = "x", ylab = "psi", dimension = 1, twodplot = persp, enhance = TRUE, 
    efactor = 0.050000000000000003, scaling.function = FALSE, type="l", ...)
{
    if(is.na(IsPowerOfTwo(resolution)))
        stop("Resolution must be a power of two")
    if(scaling.function == FALSE) {
        resolution <- resolution/2  #
#
# First obtain support widths
#
        sp <- support(filter.number = filter.number, family = family, m
             = 0, n = 0)
        lh <- c(sp$phi.lh, sp$phi.rh)
        lh <- lh[1]
        rh <- sp$phi.rh + 2 * resolution - 1
        if(verbose == TRUE)
            cat("Support of highest resolution wavelets is [", lh, 
                ", ", rh, "]\n")    #
        pic.support <- support(filter.number = filter.number, family = 
            family, m = 0, n = 0)
        pic.support <- c(pic.support$psi.lh, pic.support$psi.rh)    #
#
# Now go through all levels and see what is the lowest resolution wavelet
# that we can use to get the whole wavelet in the support range of the
# highest resolution wavelets.
#
        lowest.level <- log(resolution)/log(2)
        if(verbose == TRUE)
            cat("Lowest level is: ", lowest.level, "\n")
        selection <- NULL
        candidates <- NULL
        for(m in lowest.level:0) {
            if(verbose == TRUE) cat("Level ", m, " testing\n")  #
#
# Go through each wavelet at this level and find out
# it's support. Then check to see if it lies in the
# union of the supports of the highest resolution
# wavelets, and select it if it does.
# 
# If fact we examine all the ones that will fit, and choose one that
# is near the middle - to get a nice picture.
#
            for(n in 0:(2^(lowest.level - m) - 1)) {
                lhs <- support(filter.number = filter.number, 
                  family = family, m = m, n = n)
                rhs <- lhs$rh
                lhs <- lhs$lh
                if(verbose == TRUE)
                  cat("LHS: ", lhs, " RHS: ", rhs, "\n")
                if((lhs >= lh) && (rhs <= rh)) {
                  candidates <- c(candidates, n)
                  if(verbose == TRUE)
                    cat("Level ", m, " Position: ", n, 
                      " selected\n")
                }
            }
            if(!is.null(candidates)) {
                if(verbose == TRUE) {
                  cat("Candidates are \n")
                  print(candidates)
                }
                n <- floor(median(candidates))
                if(verbose == TRUE)
                  cat("Choosing ", n, "\n")
                selection <- list(m = m, n = n)
                lhs <- support(filter.number = filter.number, 
                  family = family, m = m, n = n)
                rhs <- lhs$rh
                lhs <- lhs$lh
                break
            }
            if(!is.null(selection))
                break
        }
#
#
#   If we haven't selected anything, then set the coefficient to
#   be one of the highest resolution coefficients. ALL of these
#   are guaranteed to be in the union of all their supports!
#   The picture will be crap though!
#
        if(is.null(selection)) selection <- list(m = 0, n = 0)  #
#
#   Build a wd object structure consisting solely of zeroes.
#
        zwd <- wd(rep(0, length = resolution * 2), filter.number = 
            filter.number, family = family, bc = "symmetric")   #
#
#   Insert a vector containing a 1 where we want to put the coefficient
#
        wd.lev <- lowest.level - selection$m
        if(verbose == TRUE)
            cat("Coefficient insertion at wd level: ", wd.lev, "\n"
                )
        if(wd.lev == 0)
            pickout <- 1
        else {
            pickout <- rep(0, 2^wd.lev)
            pickout[selection$n + 1] <- 1
        }
        zwd <- putD(zwd, level = wd.lev, v = pickout)   #
#
#   Reconstruct
#
        zwr <- wr(zwd)  #
#
#   Scales
#
        if(verbose == TRUE) {
            cat("ps: ", pic.support[1], pic.support[2], "\n")
            cat("lh,rh: ", lh, rh, "\n")
            cat("lhs,rhs: ", lhs, rhs, "\n")
        }
        aymax <- ((pic.support[2] - pic.support[1]) * (rh - lh))/(rhs - 
            lhs)
        ax <- pic.support[1] - (aymax * (lhs - lh))/(rh - lh)
        ay <- ax + aymax
        if(verbose == TRUE) cat("ax,ay ", ax, ay, "\n") #
#
#   Scale up y values, because we're actually using a higher "resolution"
#   wavelet than psi(x)
# 
        zwr <- zwr * sqrt(2)^(selection$m + 1)  #
#
#   Plot it if required
#
        x <- seq(from = ax, to = ay, length = resolution * 2)
        if(enhance == TRUE) {
            sv <- (abs(zwr) > efactor * range(abs(zwr))[2])
            sv <- (1:length(sv))[sv]
            tr <- range(sv)
            sv <- tr[1]:tr[2]
            x <- x[sv]
            zwr <- zwr[sv]
            main <- paste(main, " (Enhanced)")
        }
        if(plot.it == TRUE) {
            if(dimension == 1)
                plot(x = x, y = zwr, main = main, sub = sub, 
                  xlab = xlab, ylab = ylab, type = type, ...)
            else if(dimension == 2) {
                twodplot(x = x, y = x, z = outer(zwr, zwr), 
                  xlab = xlab, ylab = xlab, zlab = ylab, ...)
                title(main = main, sub = sub)
                invisible()
            }
            else stop("Can only do 1 or 2 dimensional plots")
        }
        else {
            if(dimension == 1)
                return(list(x = x, y = zwr))
            else if(dimension == 2)
                return(list(x = x, y = x, z = outer(zwr, zwr)))
            else stop("Can only do 1 or 2 dimensional plots")
        }
    }
    else {
        if(dimension != 1)
            stop("Can only generate one-dimensional scaling function"
                )
        if(enhance == TRUE) {
            enhance <- FALSE
            warning("Cannot enhance picture of scaling function")
        }
        if(missing(main))
            main <- "Scaling Function"
        if(missing(ylab))
            ylab <- "phi"
        if(missing(sub))
            sub <- filter.select(filter.number = filter.number, 
                family = family)$name
        phi <- ScalingFunction(filter.number = filter.number, family = 
            family, resolution = resolution)
        if(plot.it == TRUE) {
            plot(x = phi$x, y = phi$y, main = main, sub = sub, xlab
                 = xlab, ylab = ylab, type = type, ...)
        }
        else return(list(x = phi$x, y = phi$y))
    }
}
"draw.imwd"<-
function(wd, resolution = 128, ...)
{
    filter <- wd$filter
    draw.default(filter.number = filter$filter.number, family = filter$
        family, dimension = 2, resolution = resolution, ...)
}
"draw.imwdc"<-
function(wd, resolution = 128, ...)
{
    filter <- wd$filter
    draw.default(filter.number = filter$filter.number, family = filter$
        family, dimension = 2, resolution = resolution, ...)
}
"draw.mwd"<-
function(mwd, phi = 0, psi = 0, return.funct = FALSE, ...)
{
#draw.mwd
#
# plots one of the scaling or 
# wavelet functions used to create mwd
#
#check phi and psi
    if(phi > 0 && psi > 0) stop("only one of phi and psi should be nonzero"
            )
    if(phi == 0 && psi < 0)
        stop("bad psi arguement")
    if(phi < 0 && psi == 0)
        stop("bad phi arguement")
    if(phi == 0 && psi == 0)
        phi <- 1
    if(phi > mwd$filter$nphi)
        stop("There aren't that many scaling functions")
    if(psi > mwd$filter$npsi) stop("There aren't that many wavelets")   
    #for the specified case insert a single 1 and reconstruct.
    if(phi != 0) {
        main <- c("scaling function No.", phi)
        M <- matrix(rep(0, 2 * mwd$filter$nphi), nrow = mwd$filter$nphi
            )
        M[phi, 1] <- 1
        mwd$D <- matrix(rep(0, mwd$filter$npsi * mwd$fl.dbase$nvecs.d), 
            nrow = mwd$filter$npsi)
        mwd <- putC.mwd(mwd, level = 1, M)
    }
    if(psi != 0) {
        M <- matrix(rep(0, 2 * mwd$filter$npsi), nrow = mwd$filter$npsi
            )
        M[psi, 1] <- 1
        mwd$C <- matrix(rep(0, mwd$filter$nphi * mwd$fl.dbase$nvecs.c), 
            nrow = mwd$filter$nphi)
        mwd$D <- matrix(rep(0, mwd$filter$npsi * mwd$fl.dbase$nvecs.d), 
            nrow = mwd$filter$npsi)
        mwd <- putD.mwd(mwd, level = 1, M)
    }
    fun <- mwr(mwd, start.level = 1)
    x <- (2 * (0:(length(fun) - 1)))/length(fun)    #
#
#plotit
    plot(x, fun, type = "l", ...)
    if(return.funct == TRUE)
        return(fun)
}
"draw.wd"<-
function(wd, ...)
{
    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
    filter <- wd$filter
    draw.default(filter.number = filter$filter.number, family = filter$
        family, type = "l", ...)
}
"draw.wp"<-
function(wp, level, index, plot.it = TRUE, main = "Wavelet Packet", sub = paste(wp$
    name, " Level=", level, "Index= ", index), xlab = "Position", ylab = 
    "Wavelet Packet Value", ...)
{
    tmp <- drawwp.default(level = level, index = index, filter.number = wp$
        filter$filter.number, family = wp$filter$family, ...)
    if(plot.it == TRUE) {
        plot(1:length(tmp), y = tmp, main = main, sub = sub, xlab = 
            xlab, ylab = ylab, type = "l", ...)
    }
    else return(list(x = 1:length(tmp), y = tmp))
}
"draw.wst"<-
function(wst, ...)
{
    filter <- wst$filter
    draw.default(filter.number = filter$filter.number, family = filter$
        family, type = "l", ...)
}
"drawbox"<-
function(x, y, w, h, density, col)
{
    xc <- c(x, x + w, x + w, x)
    yc <- c(y, y, y + h, y + h)
    polygon(x = xc, y = yc, density = density, col = col)
}
"drawwp.default"<-
function(level, index, filter.number = 10, family = "DaubLeAsymm", resolution
     = 64 * 2^level)
{
#
# First construct a zeroed wp object 
#
    z <- rep(0, resolution) #
#
# Now take the wp transform
#
    zwp <- wp(z, filter.number = filter.number, family = family)    #
#
#
# The packet to install
#
    if(level == 0) {
        newpkt <- 1
    }
    else {
        newpkt <- rep(0, 2^level)
        newpkt[(2^level)/2] <- 1
    }
    zwp <- putpacket(zwp, level = level, index = index, packet = newpkt)    #
#
# Now set up the packet list
#
    nlev <- nlevelsWT(zwp)
    npkts <- 2^(nlev - level)
    levvec <- rep(level, npkts)
    pkt <- 0:(npkts - 1)
    basiscoef <- rep(0, npkts)
    pktlist <- list(nlevels = nlev, level = levvec, pkt = pkt)  #
#
# Do the inverse
#
    zwr <- InvBasis(zwp, pktlist = pktlist)
    zwr
}
"ewspec"<-
function(x, filter.number = 10, family = "DaubLeAsymm", UseLocalSpec = TRUE, DoSWT
     = TRUE, WPsmooth = TRUE, verbose = FALSE, smooth.filter.number = 10, 
    smooth.family = "DaubLeAsymm", smooth.levels = 3:(nlevelsWT(WPwst) - 1), 
    smooth.dev = madmad, smooth.policy = "LSuniversal", smooth.value = 0, 
    smooth.by.level = FALSE, smooth.type = "soft", smooth.verbose = FALSE, 
    smooth.cvtol = 0.01, smooth.cvnorm = l2norm, smooth.transform = I, 
    smooth.inverse = I)
{
#
#
#   Coarser is an old parameter, not needed now
#
    coarser <- 0
    if(verbose) cat("Smoothing then inversion\n")   #
#
# First compute the SWT
#
    if(DoSWT == TRUE) {
        if(verbose)
            cat("Computing nondecimated wavelet transform of data\n")
        xwdS <- wd(x, filter.number = filter.number, family = family, 
            type = "station")
    }
    else xwdS <- x
    if(UseLocalSpec == TRUE) {
        if(verbose)
            cat("Computing raw wavelet periodogram\n")
        xwdWP <- LocalSpec(xwdS, lsmooth = "none", nlsmooth = FALSE)
    }
    else xwdWP <- x
    J <- nlevelsWT(xwdWP) #
#
# Compute the vSNK matrix
#
    if(verbose)
        cat("Computing A matrix\n")
    rm <- ipndacw( - J, filter.number = filter.number, family = family) #
# Compute the inverse of the vSNK matrix
#
    if(verbose)
        cat("Computing inverse of A\n")
    irm <- solve(rm)    #
#
# Create a matrix to store the wavelet periodogram in
#
    if(verbose)
        cat("Putting wavelet periodogram into a matrix\n")
    WavPer <- matrix(0, nrow = (J - coarser), ncol = 2^J)   #
#
# Now create the Wavelet Periodogram matrix
#
#   n.b. J is coarsest  0 in wavethresh notation
#        1 is finest    J-1 in wavethresh notation
#
#   Conversion is j -> J-j
#
    for(j in 1:(J - coarser)) {
        WavPer[j,  ] <- accessD(xwdWP, lev = J - j)
    }
#
#
# Smooth the wavelet periodogram
#
    if(WPsmooth == TRUE) {
        if(verbose) {
            cat("Smoothing the wavelet periodogram\n")
            cat("Smoothing level: ")
        }
        for(j in 1:(J - coarser)) {
            if(verbose)
                cat(J - j)
            WP <- WavPer[j,  ]
            WP <- smooth.transform(WP)
            WPwst <- wst(WP, filter.number = smooth.filter.number, 
                family = smooth.family)
            if(verbose == TRUE)
                cat(".w")
            WPwstT <- threshold.wst(WPwst, levels = smooth.levels, 
                dev = smooth.dev, policy = smooth.policy, value
                 = smooth.value, by.level = smooth.by.level, 
                type = smooth.type, verbose = smooth.verbose, 
                cvtol = smooth.cvtol, cvnorm = smooth.cvnorm)
            if(verbose == TRUE)
                cat(".t")
            WPwsrR <- AvBasis(WPwstT)
            if(verbose == TRUE)
                cat(".i")
            WavPer[j,  ] <- smooth.inverse(WPwsrR)
        }
        if(verbose == TRUE)
            cat("\n")
    }
#
#
# Need a smaller inverse Rainer matrix if don't do all levels
#
    irm <- irm[1:(J - coarser), 1:(J - coarser)]    #
#
# Now multiply the inverse matrix into the WavPer
#
    S <- irm %*% WavPer #
#
# Store these levels in the xwdS object
#
    xwdS <- xwdWP
    for(j in 1:(J - coarser)) {
        xwdS <- putD(xwdS, lev = J - j, v = S[j,  ])
    }
    if(coarser > 0)
        for(j in (J - coarser + 1):J)
            xwdS <- putD(xwdS, lev = J - j, v = rep(0, 2^J))
    list(S = xwdS, WavPer = xwdWP, rm = rm, irm = irm)
}
"example.1"<-
function()
{
    x <- seq(0, 1, length = 513)
    x <- x[1:512]
    y <- rep(0, length(x))
    xsv <- (x <= 0.5)   # Left hand end
    y[xsv] <- -16 * x[xsv]^3 + 12 * x[xsv]^2
    xsv <- (x > 0.5) & (x <= 0.75)  # Middle section
    y[xsv] <- (x[xsv] * (16 * x[xsv]^2 - 40 * x[xsv] + 28))/3 - 1.5
    xsv <- x > 0.75 #Right hand end
    y[xsv] <- (x[xsv] * (16 * x[xsv]^2 - 32 * x[xsv] + 16))/3
    list(x = x, y = y)
}
"first.last"<-
function(LengthH, DataLength, type = "wavelet", bc = "periodic", current.scale
     = 0)
{
    if(type == "station" && bc != "periodic")
        stop("Can only do periodic boundary conditions with station")
    if(type != "station" && type != "wavelet")
        stop("Type can only be wavelet or station")
    levels <- log(DataLength)/log(2)
    first.last.c <- matrix(0, nrow = levels + 1, ncol = 3, dimnames = list(
        NULL, c("First", "Last", "Offset")))
    first.last.d <- matrix(0, nrow = levels - current.scale, ncol = 3, 
        dimnames = list(NULL, c("First", "Last", "Offset")))
    if(bc == "periodic") {
# Periodic boundary correction
        if(type == "wavelet") {
            first.last.c[, 1] <- rep(0, levels + 1)
            first.last.c[, 2] <- 2^(0:levels) - 1
            first.last.c[, 3] <- rev(c(0, cumsum(rev(1 + 
                first.last.c[, 2]))[1:levels]))
            first.last.d[, 1] <- rep(0, levels)
            first.last.d[, 2] <- 2^(0:(levels - 1)) - 1
            first.last.d[, 3] <- rev(c(0, cumsum(rev(1 + 
                first.last.d[, 2]))[1:(levels - 1)]))
            ntotal <- 2 * DataLength - 1
            ntotal.d <- DataLength - 1
        }
        else if(type == "station") {
            first.last.c[, 1] <- rep(0, levels + 1)
            first.last.c[, 2] <- 2^levels - 1
            first.last.c[, 3] <- rev(c(0, cumsum(rev(1 + 
                first.last.c[, 2]))[1:levels]))
            first.last.d[, 1] <- rep(0, levels)
            first.last.d[, 2] <- 2^levels - 1
            first.last.d[, 3] <- rev(c(0, cumsum(rev(1 + 
                first.last.d[, 2]))[1:(levels - 1)]))
            ntotal <- (levels + 1) * 2^levels
            ntotal.d <- levels * 2^levels
        }
    }
    else if(bc == "symmetric") {
# Symmetric boundary reflection
        first.last.c[levels + 1, 1] <- 0
        first.last.c[levels + 1, 2] <- DataLength - 1
        first.last.c[levels + 1, 3] <- 0
        ntotal <- first.last.c[levels + 1, 2] - first.last.c[levels + 1,
            1] + 1
        ntotal.d <- 0
        for(i in levels:1) {
            first.last.c[i, 1] <- trunc(0.5 * (1 - LengthH + 
                first.last.c[i + 1, 1]))
            first.last.c[i, 2] <- trunc(0.5 * first.last.c[i + 1, 2
                ])
            first.last.c[i, 3] <- first.last.c[i + 1, 3] + 
                first.last.c[i + 1, 2] - first.last.c[i + 1, 1] +
                1
            first.last.d[i, 1] <- trunc(0.5 * (first.last.c[i + 1, 
                1] - 1))
            first.last.d[i, 2] <- trunc(0.5 * (first.last.c[i + 1, 
                2] + LengthH - 2))
            if(i != levels) {
                first.last.d[i, 3] <- first.last.d[i + 1, 3] + 
                  first.last.d[i + 1, 2] - first.last.d[i + 1, 
                  1] + 1
            }
            ntotal <- ntotal + first.last.c[i, 2] - first.last.c[i, 
                1] + 1
            ntotal.d <- ntotal.d + first.last.d[i, 2] - 
                first.last.d[i, 1] + 1
        }
    }
    else if(bc == "interval") {
        first.last.d[, 1] <- rep(0, levels - current.scale)
        first.last.d[, 3] <- 2^(current.scale:(levels - 1))
        first.last.d[, 2] <- first.last.d[, 3] - 1
        first.last.c <- c(0, 2^current.scale - 1, 0)
        return(list(first.last.c = first.last.c, first.last.d = 
            first.last.d))
    }
    else {
        stop("Unknown boundary correction method")
    }
    names(ntotal) <- NULL
    names(ntotal.d) <- NULL
    list(first.last.c = first.last.c, ntotal = ntotal, first.last.d = 
        first.last.d, ntotal.d = ntotal.d)
}
"firstdot"<-
function(s)
{
    ls <- length(s)
    nc <- nchar(s)
    fd <- rep(0, ls)
    for(i in 1:ls) {
        for(j in 1:nc[i]) {
            ss <- substring(s[i], j, j)
            if(ss == ".") {
                fd[i] <- j
                break
            }
        }
    }
    fd
}
"getarrvec"<-
function(nlevels, sort = TRUE)
{
    n <- 2^nlevels
    v <- 1:n
    arrvec <- matrix(0, nrow = n, ncol = nlevels - 1)
    if(sort == TRUE) {
        for(i in 1:ncol(arrvec))
            arrvec[, i] <- sort.list(levarr(v, i))
    }
    else {
        for(i in 1:ncol(arrvec))
            arrvec[, i] <- levarr(v, i)
    }
    arrvec
}
"getpacket"<-
function(...)
UseMethod("getpacket")
"getpacket.wp"<-
function(wp, level, index, ...)
{
    if(class(wp) != "wp")
        stop("wp object is not of class wp")
    if(level > nlevelsWT(wp))
        stop("Not that many levels in wp object")
    unit <- 2^level
    LocalIndex <- unit * index + 1
    if(index > 2^(nlevelsWT(wp) - level) - 1) {
        cat("Index was too high, maximum for this level is ", 2^(wp$
            nlevels - level) - 1, "\n")
        stop("Error occured")
    }
    if(LocalIndex < 0)
        stop("Index must be  non-negative")
    packet <- wp$wp[level + 1, (LocalIndex:(LocalIndex + unit - 1))]
    packet
}
"getpacket.wpst"<-
function(wpst, level, index, ...)
{
    nlev <- nlevelsWT(wpst)
    if(level < 0)
        stop("Level must be greater than or equal to 0")
    else if(level > nlev)
        stop(paste("Level must be less than or equal to ", nlev))
    npkts <- 4^(nlev - level)
    if(index < 0)
        stop("Packet index must be greater than or equal to 0")
    else if(index > npkts - 1)
        stop(paste("Packet index must be less than or equal to ", npkts -
            1))
    pktlength <- 2^level
    lix <- 1 + wpst$avixstart[level + 1] + pktlength * index
    rix <- lix + pktlength - 1
    wpst$wpst[lix:rix]
}
"getpacket.wst"<-
function(wst, level, index, type = "D", aspect = "Identity", ...)
{
    if(type != "D" && type != "C")
        stop("Type of access must be C or D")
    class(wst) <- "wp"
    if(type == "C")
        wst$wp <- wst$Carray
    coefs <- getpacket.wp(wst, level = level, index = index)
    if(aspect == "Identity")
        return(coefs)
    else {
        fn <- get(aspect)
        return(fn(coefs))
    }
}
"getpacket.wst2D"<-
function(wst2D, level, index, type = "S", Ccode = TRUE, ...)
{
    nlev <- nlevelsWT(wst2D)
    if(level > nlev - 1)
        stop(paste("Maximum level is ", nlev - 1, " you supplied ", 
            level))
    else if(level < 0)
        stop(paste("Minimum level is 0 you supplied ", level))
    if(type != "S" && type != "H" && type != "V" && type != "D")
        stop("Type must be one of S, H, V or D")
    if(nchar(index) != nlev - level)
        stop(paste("Index must be ", nlev - level, 
            " characters long for level ", level))
    for(i in 1:nchar(index)) {
        s1 <- substring(index, i, i)
        if(s1 != "0" && s1 != "1" && s1 != "2" && s1 != "3")
            stop(paste("Character ", i, 
                " in index is not a 0, 1, 2 or 3. It is ", s1))
    }
    if(Ccode == TRUE) {
        ntype <- switch(type,
            S = 0,
            H = 1,
            V = 2,
            D = 3)
        amdim <- dim(wst2D$wst2D)
        sl <- 2^level
        out <- matrix(0, nrow = sl, ncol = sl)
        ans <- .C("getpacketwst2D",
            am = as.double(wst2D$wst2D),
            d1 = as.integer(amdim[1]),
            d12 = as.integer(amdim[1] * amdim[2]),
            maxlevel = as.integer(nlev - 1),
            level = as.integer(level),
            index = as.integer(index),
            ntype = as.integer(ntype),
            out = as.double(out),
            sl = as.integer(sl), PACKAGE = "wavethresh")
        return(matrix(ans$out, nrow = ans$sl))
    }
    else {
        x <- y <- 0
        ans <- .C("ixtoco",
            level = as.integer(level),
            maxlevel = as.integer(nlev - 1),
            index = as.integer(index),
            x = as.integer(x),
            y = as.integer(y), PACKAGE = "wavethresh")
        cellength <- 2^level
        tmpx <- switch(type,
            S = 0,
            H = 0,
            V = cellength,
            D = cellength)
        tmpy <- switch(type,
            S = 0,
            H = cellength,
            V = 0,
            D = cellength)
        x <- ans$x + tmpx + 1
        y <- ans$y + tmpy + 1
        cat("x ", x, "y: ", y, "x+cellength-1 ", x + cellength - 1, 
            "y+cellength-1", y + cellength - 1, "\n")
        return(wst2D$wst2D[level + 1, x:(x + cellength - 1), y:(y + 
            cellength - 1)])
    }
}
"guyrot"<-
function(v, n)
{
    l <- length(v)
    n <- n %% l
    if(n == 0)
        return(v)
    tmp <- v[(l - n + 1):l]
    v[(n + 1):l] <- v[1:(l - n)]
    v[1:n] <- tmp
    v
}

"image.wd"<-
function(x, strut = 10, type = "D", transform = I, ...)
{
    if(x$type != "station")
        stop("You have not supplied a nondecimated wd object")
    nlev <- nlevelsWT(x)
    if(type == "D" ) {
        m <- matrix(0, nrow = nlev, ncol = 2^nlev)
        for(i in 0:(nlev - 1)) {
            m[i,  ] <- accessD(x, lev = i)
        }
    }
    if(type == "C") {
        mC <- matrix(0, nrow = nlev + 1, ncol = 2^nlev)
        for(i in 0:nlev) {
            mC[i,  ] <- accessC(x, lev = i)
        }
    }
    nr <- nlev
    mz <- matrix(0, nrow = nlev, ncol = 2^nlev)
    if(type == "D") {
        image(transform(m[rep(1:nr, rep(strut, nr)),  ]),
            main="Wavelet coefficients")
    }
    else if(type == "C")
        image(transform(mC[rep(1:nr, rep(strut, nr)),  ]), 
             main = "Scaling function coefficients")
}
"image.wst"<-
function(x, nv, strut = 10, type = "D", transform = I, ...)
{
    m <- x$wp
    mC <- x$Carray
    nr <- nrow(m)
    nlev <- nlevelsWT(x)
    mz <- matrix(0, nrow = nrow(mC), ncol = ncol(mC))
    if(!missing(nv)) {
        pknums <- print.nv(nv, printing = FALSE)$indexlist
        mpk <- matrix(0, nrow = nrow(mC), ncol = ncol(mC))
        for(i in seq(along = pknums)) {
            lev <- nlev - i + 1
            pklength <- 2^(lev - 1)
            f <- pknums[i] * pklength + 1
            l <- f + pklength - 1
            mpk[lev, f:l] <- 1
        }
    }
    if(type == "D") {
            image(transform(m[rep(1:nr, rep(strut, nr)),  ]), 
                 main = 
                "Wavelet coefficients")
        }
    else if(type == "C")
            image(transform(mC[rep(1:nr, rep(strut, nr)),  ]), 
                 main = 
                "Scaling function coefficients"
                )
}
"imwd"<-
function(image, filter.number = 10, family = "DaubLeAsymm", type = "wavelet", 
    bc = "periodic", RetFather = TRUE, verbose = FALSE)
{
    if(verbose == TRUE)
        cat("Argument checking...")
    if(nrow(image) != ncol(image))
        stop("Number of rows and columns in image are not identical")
    if(verbose == TRUE) cat("...done\nFilter...")   #
#
#   Select wavelet filter
#
    filter <- filter.select(filter.number = filter.number, family = family)
    Csize <- nrow(image)    #
#
# Check that Csize is a power of 2
#
    nlev <- IsPowerOfTwo(Csize)
    if(is.na(nlev)) stop(paste("The image size (", Csize, 
            ") is not a power of 2"))   #
#
# Set-up first/last database
#
    if(verbose == TRUE)
        cat("...selected\nFirst/last database...")
    fl.dbase <- first.last(LengthH = length(filter$H), DataLength = Csize, 
        bc = bc, type = type)
    first.last.c <- fl.dbase$first.last.c
    first.last.d <- fl.dbase$first.last.d   #
#
# Set up answer list
#
    image.decomp <- list(nlevels = nlev, fl.dbase = fl.dbase, filter = 
        filter, type = type, bc = bc, date = date())    #
#
#
    if(verbose == TRUE) cat("...built\n")   #
#
# Ok, go round loop doing decompositions
#
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2)
    if(is.null(nbc))
        stop("Unknown boundary handling")
    if(type == "station" && bc == "symmetric")
        stop("Cannot do nondecimated transform with symmetric boundary conditions"
            )
    ntype <- switch(type,
        wavelet = 1,
        station = 2)
    if(is.null(ntype)) stop("Unknown type of transform")    #
#
#   Load up original image
#
    smoothed <- as.vector(image)
    if(verbose == TRUE) {
        cat(bc, " boundary handling\n")
        cat("Decomposing...")
    }
    for(level in seq(nrow(first.last.d), 1, -1)) {
        if(verbose == TRUE)
            cat(level - 1, "")
        LengthCin <- first.last.c[level + 1, 2] - first.last.c[level + 
            1, 1] + 1
        LengthCout <- first.last.c[level, 2] - first.last.c[level, 1] + 
            1
        LengthDout <- first.last.d[level, 2] - first.last.d[level, 1] + 
            1
        ImCC <- rep(0, (LengthCout * LengthCout))
        ImCD <- rep(0, (LengthCout * LengthDout))
        ImDC <- rep(0, (LengthDout * LengthCout))
        ImDD <- rep(0, (LengthDout * LengthDout))
        error <- 0
        z <- .C("StoIDS",
            C = as.double(smoothed),
            Csize = as.integer(LengthCin),
            firstCin = as.integer(first.last.c[level + 1, 1]),
            H = as.double(filter$H),
            LengthH = as.integer(length(filter$H)),
            LengthCout = as.integer(LengthCout),
            firstCout = as.integer(first.last.c[level, 1]),
            lastCout = as.integer(first.last.c[level, 2]),
            LengthDout = as.integer(LengthDout),
            firstDout = as.integer(first.last.d[level, 1]),
            lastDout = as.integer(first.last.d[level, 2]),
            ImCC = as.double(ImCC),
            ImCD = as.double(ImCD),
            ImDC = as.double(ImDC),
            ImDD = as.double(ImDD),
            nbc = as.integer(nbc),
            ntype = as.integer(ntype),
            error = as.integer(error), PACKAGE = "wavethresh")
        error <- z$error
        if(error != 0) {
            cat("Error was ", error, "\n")
            stop("Error reported")
        }
        smoothed <- z$ImCC
        if(RetFather == TRUE) {
            nm <- lt.to.name(level - 1, "CC")
            image.decomp[[nm]] <- z$ImCC
        }
        nm <- lt.to.name(level - 1, "CD")
        image.decomp[[nm]] <- z$ImCD
        nm <- lt.to.name(level - 1, "DC")
        image.decomp[[nm]] <- z$ImDC
        nm <- lt.to.name(level - 1, "DD")
        image.decomp[[nm]] <- z$ImDD
    }
    if(verbose == TRUE)
        cat("\nReturning answer...\n")
    image.decomp$w0Lconstant <- smoothed
    image.decomp$bc <- bc
    image.decomp$date <- date()
    class(image.decomp) <- "imwd"
    image.decomp
}
"imwr"<-
function(...)
UseMethod("imwr")
"imwr.imwd"<-
function(imwd, bc = imwd$bc, verbose = FALSE, ...)
{
    if(verbose == TRUE) cat("Argument checking...") #
#
#       Check class of imwd
#
    ctmp <- class(imwd)
    if(is.null(ctmp))
        stop("imwd has no class")
    else if(ctmp != "imwd")
        stop("imwd is not of class imwd")
    if(imwd$type == "station")
        stop("Cannot invert nonodecimated wavelet transform using imwr")
    filter <- imwd$filter
    if(verbose == TRUE)
        cat("...done\nFirst/last database...")
    fl.dbase <- imwd$fl.dbase
    first.last.c <- fl.dbase$first.last.c
    first.last.d <- fl.dbase$first.last.d
    if(verbose == TRUE)
        cat("...extracted\n")
    ImCC <- imwd$w0Lconstant
    if(verbose == TRUE) cat("Reconstructing...")    #
#
# Ok, go round loop doing reconstructions
#
    for(level in seq(2, 1 + nlevelsWT(imwd))) {
        if(verbose == TRUE)
            cat(level - 1, " ")
        LengthCin <- first.last.c[level - 1, 2] - first.last.c[level - 
            1, 1] + 1
        LengthCout <- first.last.c[level, 2] - first.last.c[level, 1] + 
            1
        LengthDin <- first.last.d[level - 1, 2] - first.last.d[level - 
            1, 1] + 1
        error <- 0
        ImOut <- rep(0, LengthCout^2)
        nbc <- switch(bc,
            periodic = 1,
            symmetric = 2)
        if(is.null(nbc))
            stop("Unknown boundary handling")
        z <- .C("StoIRS",
            ImCC = as.double(ImCC),
            ImCD = as.double(imwd[[lt.to.name(level - 2, "CD")]]),
            ImDC = as.double(imwd[[lt.to.name(level - 2, "DC")]]),
            ImDD = as.double(imwd[[lt.to.name(level - 2, "DD")]]),
            LengthCin = as.integer(LengthCin),
            firstCin = as.integer(first.last.c[level - 1, 1]),
            LengthDin = as.integer(LengthDin),
            firstDin = as.integer(first.last.d[level - 1, 1]),
            H = as.double(filter$H),
            LengthH = as.integer(length(filter$H)),
            LengthCout = as.integer(LengthCout),
            firstCout = as.integer(first.last.c[level, 1]),
            lastCout = as.integer(first.last.c[level, 2]),
            ImOut = as.double(ImOut),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
        error <- z$error
        if(error != 0) {
            cat("Error was ", error, "\n")
            stop("Error reported")
        }
# Do something with ImOut 
        ImCC <- z$ImOut
    }
    if(verbose == TRUE)
        cat("\nReturning image\n")  # Return the image
    matrix(ImCC, nrow = 2^(nlevelsWT(imwd)))
}
"imwr.imwdc"<-
function(imwd, verbose = FALSE, ...)
{
    if(verbose == TRUE)
        cat("Uncompressing...\n")
    imwd2 <- uncompress(imwd, ver = verbose)
    if(verbose == TRUE)
        cat("Reconstructing...\n")
    imwr(imwd2, verbose = verbose, ...)
}

"ipndacw"<-
function(J, filter.number = 10, family = "DaubLeAsymm", tol = 1e-100, verbose
     = FALSE, ...)
{
    if(verbose == TRUE)
        cat("Computing ipndacw\n")
    now <- proc.time()[1:2]
    if(J >= 0)
        stop("J must be negative integer")
    if(J - round(J) != 0)
        stop("J must be an integer")    #
    rmnorig <- rmname(J = J, filter.number = filter.number, family = family
        )   #
#
#   See if matrix already exists. If so, return it
#
    rm.there <- rmget(requestJ =  - J, filter.number = filter.number, 
        family = family)
    if(!is.null(rm.there)) {
        if(verbose == TRUE)
            cat("Returning precomputed version: using ", rm.there, 
                "\n")
        speed <- proc.time()[1:2] - now
        if(verbose == TRUE)
            cat("Took ", sum(speed), " seconds\n")
        rmnexists <- rmname(J =  - rm.there, filter.number = 
            filter.number, family = family)
        tmp <- get(rmnexists, envir=WTEnv)[1:( - J), 1:( - J)]
        assign(rmnorig, tmp, envir=WTEnv)
        return(tmp)
    }
#
#
#   See if partially computed matrix exists. If so, use it.
#
    if(J != -1) {
        for(j in (1 + J):(-1)) {
            rmn <- rmname(J = j, filter.number = filter.number, 
                family = family)
            if(exists(rmn, envir=WTEnv)) {
                if(verbose == TRUE) {
                  cat("Partial matrix: ", rmn, " exists (")
                  cat(paste(round(100 - (100 * (j * j))/(J * J),
                    digits = 1), "% left to do)\n", sep = ""))
                }
                fmat <- rep(0, J * J)
                H <- filter.select(filter.number = 
                  filter.number, family = family)$H
                error <- 0
                answer <- .C("rainmatPARTIAL",
                  J = as.integer( - J),
                  j = as.integer( - j),
                  H = as.double(H),
                  LengthH = as.integer(length(H)),
                  fmat = as.double(fmat),
                  tol = as.double(tol),
                  error = as.integer(error), PACKAGE = "wavethresh")
                if(answer$error != 0)
                  stop(paste("Error code was ", answer$error))
                m <- matrix(answer$fmat, nrow =  - J)
                m[1:( - j), 1:( - j)] <- get(rmn, envir=WTEnv)
                nm <- as.character(-1:J)
                dimnames(m) <- list(nm, nm)
                speed <- proc.time()[1:2] - now
                if(verbose == TRUE)
                  cat("Took ", sum(speed), " seconds\n")
                assign(rmnorig, m, envir=WTEnv)
                return(m)
            }
        }
    }
#
#
#   Otherwise have to compute whole matrix
#
    fmat <- rep(0, J * J)
    H <- filter.select(filter.number = filter.number, family = family)$H
    error <- 0
    answer <- .C("rainmatPARENT",
        J = as.integer( - J),
        H = as.double(H),
        LengthH = as.integer(length(H)),
        fmat = as.double(fmat),
        tol = as.double(tol),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(answer$error != 0)
        stop(paste("Error code was ", answer$error))
    speed <- proc.time()[1:2] - now
    if(verbose == TRUE)
        cat("Took ", sum(speed), " seconds\n")
    m <- matrix(answer$fmat, nrow =  - J)
    nm <- as.character(-1:J)
    dimnames(m) <- list(nm, nm)
    assign(rmnorig, m, envir=WTEnv)
    m
}
"irregwd"<-
function(gd, filter.number = 2, family = "DaubExPhase", bc = "periodic", 
    verbose = FALSE)
{
    type <- "wavelet"
    if(verbose == TRUE)
        cat("wd: Argument checking...")
    ctmp <- class(gd)
    if(is.null(ctmp))
        stop("gd has no class")
    else if(ctmp != "griddata")
        stop("gd is not of class griddata")
    data <- gd$gridy
    if(!is.atomic(data))
        stop("Data is not atomic")
    DataLength <- length(data)  #
#
# Check that we have a power of 2 data elements
#
    nlevels <- nlevelsWT(data)    #
    if(is.na(nlevels)) stop("Data length is not power of two")  
    # Check for correct type
#
    if(type != "wavelet" && type != "station")
        stop("Unknown type of wavelet decomposition")
    if(type == "station" && bc != "periodic") stop(
            "Can only do periodic boundary conditions with station"
            )   #
# Select the appropriate filter
#
    if(verbose == TRUE)
        cat("...done\nFilter...")
    filter <- filter.select(filter.number = filter.number, family = family)
        #
#
# Build the first/last database
#
    if(verbose == TRUE)
        cat("...selected\nFirst/last database...")
    fl.dbase <- first.last(LengthH = length(filter$H), DataLength = 
        DataLength, type = type, bc = bc)   #
#
# Save time series attribute if there is one
#
    dtsp <- tsp(data)   #
#
# Put in the data
#
    C <- rep(0, fl.dbase$ntotal)
    C[1:DataLength] <- data #
    if(verbose == TRUE)
        error <- 1
    else error <- 0
    if(verbose == TRUE) cat("built\n")  #
#
# Compute the decomposition
#
    if(verbose == TRUE)
        cat("Decomposing...\n")
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2)
    if(is.null(nbc))
        stop("Unknown boundary condition")
    ntype <- switch(type,
        wavelet = 1,
        station = 2)
    if(is.null(filter$G)) {
        wavelet.decomposition <- .C("wavedecomp",
            C = as.double(C),
            D = as.double(rep(0, fl.dbase$ntotal.d)),
            H = as.double(filter$H),
            LengthH = as.integer(length(filter$H)),
            nlevels = as.integer(nlevels),
            firstC = as.integer(fl.dbase$first.last.c[, 1]),
            lastC = as.integer(fl.dbase$first.last.c[, 2]),
            offsetC = as.integer(fl.dbase$first.last.c[, 3]),
            firstD = as.integer(fl.dbase$first.last.d[, 1]),
            lastD = as.integer(fl.dbase$first.last.d[, 2]),
            offsetD = as.integer(fl.dbase$first.last.d[, 3]),
            ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
        tmp <- .C("computec",
            n = as.integer(length(gd$Gleft)),
            c = as.double(rep(0, fl.dbase$ntotal.d)),
            gridn = as.integer(length(gd$G)),
            G = as.double(gd$G),
            Gindex = as.integer(gd$Gindex),
            H = as.double(filter$H),
            LengthH = as.integer(length(filter$H)),
            nbc = as.integer(nbc), PACKAGE = "wavethresh")
    }
    else {
        wavelet.decomposition <- .C("comwd",
            CR = as.double(Re(C)),
            CI = as.double(Im(C)),
            LengthC = as.integer(fl.dbase$ntotal),
            DR = as.double(rep(0, fl.dbase$ntotal.d)),
            DI = as.double(rep(0, fl.dbase$ntotal.d)),
            LengthD = as.integer(fl.dbase$ntotal.d),
            HR = as.double(Re(filter$H)),
            HI = as.double( - Im(filter$H)),
            GR = as.double(Re(filter$G)),
            GI = as.double( - Im(filter$G)),
            LengthH = as.integer(length(filter$H)),
            nlevels = as.integer(nlevels),
            firstC = as.integer(fl.dbase$first.last.c[, 1]),
            lastC = as.integer(fl.dbase$first.last.c[, 2]),
            offsetC = as.integer(fl.dbase$first.last.c[, 3]),
            firstD = as.integer(fl.dbase$first.last.d[, 1]),
            lastD = as.integer(fl.dbase$first.last.d[, 2]),
            offsetD = as.integer(fl.dbase$first.last.d[, 3]),
            ntype = as.integer(ntype),
            nbc = as.integer(nbc),
            error = as.integer(error), PACKAGE = "wavethresh")
    }
    if(verbose == TRUE)
        cat("done\n")
    error <- wavelet.decomposition$error
    if(error != 0) {
        cat("Error ", error, " occured in wavedecomp\n")
        stop("Error")
    }
    if(is.null(filter$G)) {
        l <- list(C = wavelet.decomposition$C, D = 
            wavelet.decomposition$D, c = tmp$c * (tmp$c > 0), 
            nlevels = nlevelsWT(wavelet.decomposition), fl.dbase = 
            fl.dbase, filter = filter, type = type, bc = bc, date
             = date())
    }
    else {
        l <- list(C = complex(real = wavelet.decomposition$CR,
		imaginary = 
            wavelet.decomposition$CI), D = complex(real = 
            wavelet.decomposition$DR, imaginary = wavelet.decomposition$DI
            ), nlevels = nlevelsWT(wavelet.decomposition), fl.dbase = 
            fl.dbase, filter = filter, type = type, bc = bc, date
             = date())
    }
    class(l) <- "irregwd"
    if(!is.null(dtsp))
        tsp(l) <- dtsp
    l
}
"l2norm"<-
function(u, v)
sqrt(sum((u - v)^2))

"levarr"<-
function(v, levstodo)
{
    if(levstodo != 0) {
        sv <- seq(from = 1, to = length(v), by = 2)
        return(c(levarr(v[sv], levstodo - 1), levarr(v[ - sv], levstodo -
            1)))
    }
    else return(v)
}
"linfnorm"<-
function(u, v)
{
    max(abs(u - v))
}
"lt.to.name"<-
function(level, type)
{
#
# This function converts the level and type (horizontal, vertical, diagonal)
# of wavelet coefficients to a character string "wnLx" which should be
# interpreted as "nth Level, coefficients x", where x is 1, 2 or 3 in the
# scheme of Mallat. (So 1 is horizontal, 2 is vertical and 3 is diagonal).
# w is on the front to indicate that these are wavelet coefficients
#
    return(paste("w", as.character(level), "L", switch(type,
        CD = "1",
        DC = "2",
        DD = "3",
        CC = "4"), sep = ""))
}
"madmad"<-
function(x)
mad(x)^2
"makegrid"<-
function(t, y, gridn = 2^(floor(log(length(t) - 1, 2)) + 1))
{
    tmp <- .C("makegrid",
        x = as.double(t),
        y = as.double(y),
        n = length(t),
        gridt = as.double(rep(0, gridn)),
        gridy = as.double(rep(0, gridn)),
        gridn = as.integer(gridn),
        G = as.double(rep(0, gridn)),
        Gindex = as.integer(rep(0, gridn)), PACKAGE = "wavethresh")
    l <- list(gridt = tmp$gridt, gridy = tmp$gridy, G = tmp$G, Gindex = tmp$
        Gindex)
    class(l) <- "griddata"
    l
}
"makewpstDO"<-
function(timeseries, groups, filter.number = 10, family = "DaubExPhase", mincor
     = 0.69999999999999996)
{
#
#
# Using the data in time series (which should be a length a power of two)
# and the group information (only two groups presently). Create an object
# of class wpstDO (nondecimated wavelet packet Discrimination Object).
#
# Given this wpstDO and another timeseries a function exists to predict
# the group membership of each timeseries element
#
#
# First build nondecimated wavelet packet object
#
    twpst <- wpst(timeseries, filter.number = filter.number, family = 
        family) #
#
# Now convert this to a w2d object including the group information.
#
    tw2d <- wpst2discr(wpstobj = twpst, groups = groups)   #
#
# Now extract the best 1D classifying columns.
#
    tBP <- Best1DCols(w2d = tw2d, mincor = mincor)  #
#
# Do a discriminant analysis
#
    tBPd <- BMdiscr(tBP)
    l <- list(BPd = tBPd, BP = tBP, filter = twpst$filter)
    class(l) <- "wpstDO"
    l
}
"mfilter.select"<-
function(type = "Geronimo")
{
#
# mfilter.select
# returns the filter information for a specified
# multiple wavelet basis
#
# Copyright Tim Downie 1995-6.
#
#
    if(type == "Geronimo") {
        name <- "Geronimo Multiwavelets"
        nphi <- 2
        npsi <- 2
        NH <- 4
        ndecim <- 2
        H <- rep(0, 16)
        G <- rep(0, 16)
        H[1] <- 0.42426406871193001
        H[2] <- 0.80000000000000004
        H[3] <- -0.050000000000000003
        H[4] <- -0.21213203435596001
        H[5] <- 0.42426406871193001
        H[7] <- 0.45000000000000001
        H[8] <- 0.70710678118655002
        H[11] <- 0.45000000000000001
        H[12] <- -0.21213203435596001
        H[15] <- -0.050000000000000003  #
# H6,9,10,13,14,16 are zero.
#
        G[1] <- -0.050000000000000003
        G[2] <- -0.21213203435596401
        G[3] <- 0.070710678118654793
        G[4] <- 0.29999999999999999
        G[5] <- 0.45000000000000001
        G[6] <- -0.70710678118654802
        G[7] <- -0.63639610306789296
        G[9] <- 0.45000000000000001
        G[10] <- -0.21213203435596401
        G[11] <- 0.63639610306789296
        G[12] <- -0.29999999999999999
        G[13] <- -0.050000000000000003
        G[15] <- -0.070710678118654793  #
# G8,14,16 are zero.
#
    }
    else if(type == "Donovan3") {
        name <- "Donovan Multiwavelets, 3 functions"
        nphi <- 3
        npsi <- 3
        NH <- 4
        ndecim <- 2
        H <- rep(0, 36)
        G <- rep(0, 36)
        H[2] <- ( - sqrt(154) * (3 + 2 * sqrt(5)))/3696
        H[3] <- (sqrt(14) * (2 + 5 * sqrt(5)))/1232
        H[10] <- ( - sqrt(2) * (3 + 2 * sqrt(5)))/44
        H[11] <- (sqrt(154) * (67 + 30 * sqrt(5)))/3696
        H[12] <- (sqrt(14) * (-10 + sqrt(5)))/112
        H[19] <- 1/sqrt(2)
        H[20] <- (sqrt(154) * (67 - 30 * sqrt(5)))/3696
        H[21] <- (sqrt(14) * (10 + sqrt(5)))/112
        H[23] <- (3 * sqrt(2))/8
        H[24] <- (sqrt(22) * (-4 + sqrt(5)))/88
        H[26] <- (sqrt(22) * (32 + 7 * sqrt(5)))/264
        H[27] <- (sqrt(2) * (-5 + 4 * sqrt(5)))/88
        H[28] <- (sqrt(2) * (-3 + 2 * sqrt(5)))/44
        H[29] <- (sqrt(154) * (-3 + 2 * sqrt(5)))/3696
        H[30] <- (sqrt(14) * (-2 + 5 * sqrt(5)))/1232
        H[31] <- sqrt(154)/22
        H[32] <- (3 * sqrt(2))/8
        H[33] <- (sqrt(22) * (4 + sqrt(5)))/88
        H[34] <-  - sqrt(70)/22
        H[35] <- (sqrt(22) * (-32 + 7 * sqrt(5)))/264
        H[36] <- ( - sqrt(2) * (5 + 4 * sqrt(5)))/88    #
# H1,4,5,6,7,8,9,13,14,15,16,17,18,22,25 are zero.
#
        G[5] <- (sqrt(154) * (3 + 2 * sqrt(5)))/3696
        G[6] <- ( - sqrt(14) * (2 + 5 * sqrt(5)))/1232
        G[8] <- ( - sqrt(7) * (1 + sqrt(5)))/336
        G[9] <- (sqrt(77) * (-1 + 3 * sqrt(5)))/1232
        G[13] <- (sqrt(2) * (3 + 2 * sqrt(5)))/44
        G[14] <- ( - sqrt(154) * (67 + 30 * sqrt(5)))/3696
        G[15] <- (sqrt(14) * (10 - sqrt(5)))/112
        G[16] <- ( - sqrt(11) * (1 + sqrt(5)))/44
        G[17] <- (sqrt(7) * (29 + 13 * sqrt(5)))/336
        G[18] <- (sqrt(77) * (-75 + 17 * sqrt(5)))/1232
        G[20] <- (sqrt(77) * (-2 + sqrt(5)))/264
        G[21] <- (sqrt(7) * (13 - 6 * sqrt(5)))/88
        G[22] <- 1/sqrt(2)
        G[23] <- (sqrt(154) * (-67 + 30 * sqrt(5)))/3696
        G[24] <- ( - sqrt(14) * (10 + sqrt(5)))/112
        G[26] <- (sqrt(7) * (-29 + 13 * sqrt(5)))/336
        G[27] <- ( - sqrt(77) * (75 + 17 * sqrt(5)))/1232
        G[28] <- 13/22
        G[29] <- ( - sqrt(77) * (2 + sqrt(5)))/264
        G[30] <- ( - sqrt(7) * (13 + 6 * sqrt(5)))/88
        G[31] <- (sqrt(2) * (3 - 2 * sqrt(5)))/44
        G[32] <- (sqrt(154) * (3 - 2 * sqrt(5)))/3696
        G[33] <- (sqrt(14) * (2 - 5 * sqrt(5)))/1232
        G[34] <- (sqrt(11) * (1 - sqrt(5)))/44
        G[35] <- (sqrt(7) * (1 - sqrt(5)))/336
        G[36] <- ( - sqrt(77) * (1 + 3 * sqrt(5)))/1232 #
# G1,2,3,4,7,10,11,12,19,25 are zero.
#
    }
    else (stop("bad filter specified\n"))
    return(list(type = type, name = name, nphi = nphi, npsi = npsi, NH = NH,
        ndecim = ndecim, H = H, G = G))
}
"mfirst.last"<-
function(LengthH, nlevels, ndecim, type = "wavelet", bc = "periodic")
{
#
# mfirst.last
# Sets up a coefficient data base for a multiple wavelet object
# The structure is analogous to that used in first.last
# but returns more information required by mwd and mwr.
#
# Copyright  Tim Downie 1995-1996
#
# 
    if(type != "wavelet") stop("Type can only be wavelet")
    first.last.c <- matrix(0, nrow = nlevels + 1, ncol = 3, dimnames = list(
        NULL, c("First", "Last", "Offset")))
    first.last.d <- matrix(0, nrow = nlevels, ncol = 3, dimnames = list(
        NULL, c("First", "Last", "Offset")))
    if(bc == "periodic") {
# Periodic boundary correction
        if(type == "wavelet") {
            first.last.c[, 1] <- rep(0, nlevels + 1)
            first.last.c[, 2] <- ndecim^(0:nlevels) - 1
            first.last.c[, 3] <- rev(c(0, cumsum(rev(1 + 
                first.last.c[, 2]))[1:nlevels]))
            first.last.d[, 1] <- rep(0, nlevels)
            first.last.d[, 2] <- ndecim^(0:(nlevels - 1)) - 1
            first.last.d[, 3] <- rev(c(0, cumsum(rev(1 + 
                first.last.d[, 2]))[1:(nlevels - 1)]))
            nvecs.c <- first.last.c[1, 3] + 1
            nvecs.d <- first.last.d[1, 3] + 1
        }
        else if(type == "station") {
#
#
# in case nondecimated Multiple wavelet transform is implemented
# then this code might be of use (will need adapting)
# 
            first.last.c[, 1] <- rep(0, nlevels + 1)
            first.last.c[, 2] <- 2^nlevels - 1
            first.last.c[, 3] <- rev(c(0, cumsum(rev(1 + 
                first.last.c[, 2]))[1:nlevels]))
            first.last.d[, 1] <- rep(0, nlevels)
            first.last.d[, 2] <- 2^nlevels - 1
            first.last.d[, 3] <- rev(c(0, cumsum(rev(1 + 
                first.last.d[, 2]))[1:(nlevels - 1)]))
            ntotal <- (nlevels + 1) * 2^nlevels
            ntotal.d <- nlevels * 2^nlevels
        }
    }
    else if(bc == "symmetric") {
# Symmetric boundary reflection
        first.last.c[nlevels + 1, 1] <- 0
        first.last.c[nlevels + 1, 2] <- 2^nlevels - 1
        first.last.c[nlevels + 1, 3] <- 0
        nvecs.c <- first.last.c[nlevels + 1, 2] - first.last.c[nlevels + 
            1, 1] + 1
        nvecs.d <- 0
        for(i in nlevels:1) {
            first.last.c[i, 1] <- trunc(0.5 * (1 - LengthH + 
                first.last.c[i + 1, 1]))
            first.last.c[i, 2] <- trunc(0.5 * first.last.c[i + 1, 2
                ])
            first.last.c[i, 3] <- first.last.c[i + 1, 3] + 
                first.last.c[i + 1, 2] - first.last.c[i + 1, 1] +
                1
            first.last.d[i, 1] <- trunc(0.5 * (first.last.c[i + 1, 
                1] - 1))
            first.last.d[i, 2] <- trunc(0.5 * (first.last.c[i + 1, 
                2] + LengthH - 2))
            if(i != nlevels) {
                first.last.d[i, 3] <- first.last.d[i + 1, 3] + 
                  first.last.d[i + 1, 2] - first.last.d[i + 1, 
                  1] + 1
            }
            nvecs.c <- nvecs.c + first.last.c[i, 2] - first.last.c[
                i, 1] + 1
            nvecs.d <- nvecs.d + first.last.d[i, 2] - first.last.d[
                i, 1] + 1
        }
    }
    else {
        stop("Unknown boundary correction method")
    }
    names(nvecs.c) <- NULL
    names(nvecs.d) <- NULL
    list(first.last.c = first.last.c, nvecs.c = nvecs.c, first.last.d = 
        first.last.d, nvecs.d = nvecs.d)
}
"modernise"<-
function(...)
UseMethod("modernise")
"modernise.wd"<-
function(wd, ...)
{
    if(IsEarly(wd)) {
        cat("Converting wavelet object to latest release\n")
        wd$type <- "wavelet"
        wd$date <- date()
    }
    else cat("Object is already up to date\n")
    wd
}
"mpostfilter"<-
function(C, prefilter.type, filter.type, nphi, npsi, ndecim, nlevels, verbose
     = FALSE)
{
    ndata <- ndecim^nlevels * nphi
    if(prefilter.type == "Repeat")
        ndata <- ndecim^(nlevels - 1) * nphi
    data <- rep(0, ndata)
    if(filter.type == "Geronimo") {
        if(prefilter.type == "Minimal") {
            if(verbose == TRUE)
                cat(" O.K.\nPostfilter (Minimal)\n")
            w <- 1
            data[(1:(ndata/2)) * 2 - 1] <- 2/w * C[2, (1:(ndata/2))
                ]
            data[(1:(ndata/2)) * 2] <-  - sqrt(2)/w * C[1, (1:(
                ndata/2))] + 4/w * C[2, (1:(ndata/2))]
        }
        else if(prefilter.type == "Identity") {
            if(verbose == TRUE)
                cat(" O.K.\nPostfilter (identity)\n")
            data[(1:(ndata/2)) * 2 - 1] <- C[1, (1:(ndata/2))]
            data[(1:(ndata/2)) * 2] <- C[2, (1:(ndata/2))]
        }
        else if(prefilter.type == "Repeat") {
            if(verbose == TRUE)
                cat(" O.K.\nPostfilter (weighted average)\n")
            for(k in 1:ndata)
                data[k] <- (C[2, k] + C[1, k]/sqrt(2))/2
        }
        else if(prefilter.type == "Interp" || prefilter.type == 
            "default") {
            if(verbose == TRUE)
                cat(" O.K.\nPostfilter (interpolation)\n")
            t <- sqrt(96/25)
            u <- sqrt(3)
            data[2 * (1:(ndata/2))] <- u * C[2, (1:(ndata/2))]
            data[2 * (2:(ndata/2)) - 1] <- t * C[1, (2:(ndata/2))] - 
                0.29999999999999999 * (data[2 * (2:(ndata/2)) - 
                2] + data[2 * (2:(ndata/2))])
            data[1] <- t * C[1, 1] - 0.29999999999999999 * (data[
                ndata] + data[2])
        }
        else if(prefilter.type == "Xia") {
            if(verbose == TRUE)
                cat(" O.K.\nPostfilter (Xia)\n")
            epsilon1 <- 0
            epsilon2 <- 0.10000000000000001
            root2 <- sqrt(2)
            x <- (2 * root2)/(5 * (root2 * epsilon2 - epsilon1))
            a <- (x - epsilon1 + epsilon2 * 2 * root2)/2
            b <- (x + epsilon1 - epsilon2 * 2 * root2)/2
            c <- (x + 4 * epsilon1 - epsilon2 * 3 * root2)/(root2 * 
                2)
            d <- (x - 4 * epsilon1 + epsilon2 * 3 * root2)/(root2 * 
                2)
            data[2 * (1:(ndata/2))] <- d * C[1, 1:(ndata/2)] - b * 
                C[2, 1:(ndata/2)]
            data[2 * (1:(ndata/2)) - 1] <- a * C[2, 1:(ndata/2)] - 
                c * C[1, 1:(ndata/2)]
        }
        else if(prefilter.type == "Roach1") {
            q1 <- 0.32982054290000001
            q2 <- 0.23184851840000001
            q3 <- 0.8187567536
            q4 <- -0.29459505809999997
            q5 <- -0.1629787369
            q6 <- 0.23184851840000001
            q7 <- -0.23184851840000001
            q8 <- -0.1629787369
            q9 <- 0.29459505809999997
            q10 <- 0.8187567536
            q11 <- -0.23184851840000001
            q12 <- 0.32982054290000001
            nn <- (ndata - 2)/2
            QB <- matrix(c(q2, q1, q8, q7), ncol = 2, byrow = TRUE)
            QA <- matrix(c(q4, q3, q10, q9), ncol = 2, byrow = TRUE)
            QZ <- matrix(c(q6, q5, q12, q11), ncol = 2, byrow = TRUE)
            partition <- matrix(data, nrow = 2, byrow = FALSE)
            partition[, (2:nn)] <- QB %*% C[, (2:nn) - 1] + QA %*% 
                C[, (2:nn)] + QZ %*% C[, (2:nn) + 1]
            partition[, 1] <- QB %*% C[, nn + 1] + QA %*% C[, 1] + 
                QZ %*% C[, 2]
            partition[, nn + 1] <- QB %*% C[, nn] + QA %*% C[, nn + 
                1] + QZ %*% C[, 1]
            data <- c(partition)
        }
        else if(prefilter.type == "Roach3") {
            q1 <- 0.084397403440000004
            q2 <- -0.0036003129089999999
            q3 <- 0.084858161210000005
            q4 <- 0.99279918550000001
            q5 <- -0.00015358592229999999
            q6 <- -0.0036003129089999999
            q7 <- -0.0036003129089999999
            q8 <- 0.00015358592229999999
            q9 <- 0.99279918550000001
            q10 <- -0.084858161210000005
            q11 <- -0.0036003129089999999
            q12 <- -0.084397403440000004
            nn <- (ndata - 2)/2
            QZ <- matrix(c(q7, q8, q1, q2), ncol = 2, byrow = TRUE)
            QA <- matrix(c(q9, q10, q3, q4), ncol = 2, byrow = TRUE)
            QB <- matrix(c(q11, q12, q5, q6), ncol = 2, byrow = TRUE)
            partition <- matrix(data, nrow = 2, byrow = FALSE)
            partition[, (2:nn)] <- QB %*% C[, (2:nn) - 1] + QA %*% 
                C[, (2:nn)] + QZ %*% C[, (2:nn) + 1]
            partition[, 1] <- QB %*% C[, nn + 1] + QA %*% C[, 1] + 
                QZ %*% C[, 2]
            partition[, nn + 1] <- QB %*% C[, nn] + QA %*% C[, nn + 
                1] + QZ %*% C[, 1]
            data <- c(partition)
        }
        else stop("Specified postfilter not available for given multiwavelet"
                )
    }
    else if(filter.type == "Donovan3") {
        if(prefilter.type == "Identity") {
            if(verbose == TRUE)
                cat(" O.K.\nPostfilter (identity)\n")
            data[(1:(ndata/3)) * 3 - 2] <- C[1, (1:(ndata/3))]
            data[(1:(ndata/3)) * 3 - 1] <- C[2, (1:(ndata/3))]
            data[(1:(ndata/3)) * 3] <- C[3, (1:(ndata/3))]
        }
        else if(prefilter.type == "Linear") {
            cat(" O.K.\nPostfilter (Linear)\n")
            if(verbose == TRUE)
                data[(1:(ndata/3)) * 3 - 2] <- C[1, (1:(ndata/3
                  ))] * -4.914288 + 4.914288 * C[2, (1:(ndata/3
                  ))]
            data[(1:(ndata/3)) * 3 - 1] <- C[1, (1:(ndata/3))] * 
                -2.778375 + 3.778375 * C[2, (1:(ndata/3))]
            data[(1:(ndata/3)) * 3] <- C[1, (1:(ndata/3))] * 
                -2.298365 + 3.298365 * C[2, (1:(ndata/3))] + C[
                3, (1:(ndata/3))]
        }
        else if(prefilter.type == "Interp" || prefilter.type == 
            "default") {
            if(verbose == TRUE)
                cat(" O.K.\nPostfilter (interpolation)\n")
            w <- sqrt(5)
            lc <- length(data)/3
            data[3 * (0:(lc - 1)) + 1] <- C[1, 1:lc] * sqrt(11/7)
            data[2] <- ( - (2 + 6 * w) * C[1, lc] - (3 + 2 * w) * C[
                1, 1] + 6 * sqrt(77) * C[2, 1] + ((103 - 24 * w
                ) * sqrt(7))/(16 - 5 * w) * C[3, 1])/(9 * sqrt(
                77))
            data[3 * (1:(lc - 1)) + 2] <- ( - (2 + 6 * w) * C[1, 1:(
                lc - 1)] - (3 + 2 * w) * C[1, (2:lc)] + 6 * 
                sqrt(77) * C[2, (2:lc)] + ((103 - 24 * w) * 
                sqrt(7))/(16 - 5 * w) * C[3, (2:lc)])/(9 * sqrt(
                77))
            data[3] <- ((-3 + 2 * w)/(3 * sqrt(231)) * C[1, lc] + (
                -2
                 + 6 * w)/(3 * sqrt(231)) * C[1, 1] + 2/sqrt(3) *
                C[2, 1] + (306 - 112 * w)/((16 - 5 * w) * 3 * 
                sqrt(33)) * C[3, 1])/sqrt(3)
            data[3 * (2:lc)] <- ((-3 + 2 * w)/(3 * sqrt(231)) * C[1,
                (1:(lc - 1))] + (-2 + 6 * w)/(3 * sqrt(231)) * 
                C[1, (2:lc)] + 2/sqrt(3) * C[2, (2:lc)] + (306 - 
                112 * w)/((16 - 5 * w) * 3 * sqrt(33)) * C[3, (
                2:lc)])/sqrt(3)
        }
        else stop("Specified postfilter not available for given multiwavelet"
                )
    }
    else stop("No postfilters for type of multiwavelet")
    return(data)
}
"mprefilter"<-
function(data, prefilter.type, filter.type, nlevels, nvecs.c, nphi, npsi, 
    ndecim, verbose = FALSE)
{
#function that takes original data and computes the starting level
#coefficients for the wavelet decompostion
#
    ndata <- length(data)
    C <- matrix(rep(0, nvecs.c * nphi), nrow = nphi)    #
#jump to type of multiwavelet
    if(filter.type == "Geronimo") {
        if(prefilter.type == "Minimal") {
            if(verbose == TRUE)
                cat("  O.K.\nPrefilter (Minimal)...")
            w <- 1
            C[1, 1:(ndata/2)] <- w * sqrt(2) * data[(1:(ndata/2)) * 
                2 - 1] - w/sqrt(2) * data[(1:(ndata/2)) * 2]
            C[2, 1:(ndata/2)] <- w * 0.5 * data[(1:(ndata/2)) * 2 - 
                1]
        }
        else if(prefilter.type == "Identity") {
            if(verbose == TRUE)
                cat("  O.K.\nPrefilter (Identity)...")
            for(l in 1:nphi) {
                C[l, 1:(ndata/nphi)] <- data[(0:((ndata/nphi) - 
                  1)) * nphi + l]
            }
        }
        else if(prefilter.type == "Repeat") {
            if(verbose == TRUE)
                cat("  O.K.\nRepeating signal...")
            C[1, 1:(ndata)] <- data[1:ndata] * sqrt(2)
            C[2, 1:(ndata)] <- data[1:ndata]
        }
        else if(prefilter.type == "Interp" || prefilter.type == 
            "default") {
            if(verbose == TRUE)
                cat("  O.K.\nPrefilter (interpolation)...")
            r <- sqrt(25/96)
            s <- sqrt(1/3)
            a <- -0.29999999999999999
            C[2, (1:(ndata/2))] <- s * data[2 * (1:(ndata/2))]
            C[1, 1] <- r * (data[1] - a * (data[ndata] + data[2]))
            C[1, (2:(ndata/2))] <- r * (data[2 * (2:(ndata/2)) - 1] -
                a * (data[2 * (2:(ndata/2)) - 2] + data[2 * (2:(
                ndata/2))]))
        }
        else if(prefilter.type == "Xia") {
            if(verbose == TRUE)
                cat("  O.K.\nPrefilter (Xia) ...")
            epsilon1 <- 0
            epsilon2 <- 0.10000000000000001
            root2 <- sqrt(2)
            x <- (2 * root2)/(5 * (root2 * epsilon2 - epsilon1))
            a <- (x - epsilon1 + epsilon2 * 2 * root2)/2
            b <- (x + epsilon1 - epsilon2 * 2 * root2)/2
            c <- (x + 4 * epsilon1 - epsilon2 * 3 * root2)/(root2 * 
                2)
            d <- (x - 4 * epsilon1 + epsilon2 * 3 * root2)/(root2 * 
                2)
            C[1, (1:(ndata/2))] <- a * data[2 * (1:(ndata/2))] + b * 
                data[2 * (1:(ndata/2)) - 1]
            C[2, (1:(ndata/2))] <- c * data[2 * (1:(ndata/2))] + d * 
                data[2 * (1:(ndata/2)) - 1]
        }
        else if(prefilter.type == "Roach1") {
            q1 <- 0.32982054290000001
            q2 <- 0.23184851840000001
            q3 <- 0.8187567536
            q4 <- -0.29459505809999997
            q5 <- -0.1629787369
            q6 <- 0.23184851840000001
            q7 <- -0.23184851840000001
            q8 <- -0.1629787369
            q9 <- 0.29459505809999997
            q10 <- 0.8187567536
            q11 <- -0.23184851840000001
            q12 <- 0.32982054290000001
            QB <- matrix(c(q2, q1, q8, q7), ncol = 2, byrow = TRUE)
            QA <- matrix(c(q4, q3, q10, q9), ncol = 2, byrow = TRUE)
            QZ <- matrix(c(q6, q5, q12, q11), ncol = 2, byrow = TRUE)
            nn <- (ndata - 2)/2
            partition <- matrix(data, nrow = 2, byrow = FALSE)
            C[, (2:nn)] <- QB %*% partition[, (2:nn) - 1] + QA %*% 
                partition[, (2:nn)] + QZ %*% partition[, (2:nn) +
                1]
            C[, 1] <- QB %*% partition[, nn + 1] + QA %*% partition[
                , 1] + QZ %*% partition[, 2]
            C[, nn + 1] <- QB %*% partition[, nn] + QA %*% 
                partition[, nn + 1] + QZ %*% partition[, 1]
        }
        else if(prefilter.type == "Roach3") {
            q1 <- 0.084397403440000004
            q2 <- -0.0036003129089999999
            q3 <- 0.084858161210000005
            q4 <- 0.99279918550000001
            q5 <- -0.00015358592229999999
            q6 <- -0.0036003129089999999
            q7 <- -0.0036003129089999999
            q8 <- 0.00015358592229999999
            q9 <- 0.99279918550000001
            q10 <- -0.084858161210000005
            q11 <- -0.0036003129089999999
            q12 <- -0.084397403440000004
            nn <- (ndata - 2)/2
            QB <- matrix(c(q7, q8, q1, q2), ncol = 2, byrow = FALSE)
            QA <- matrix(c(q9, q10, q3, q4), ncol = 2, byrow = FALSE)
            QZ <- matrix(c(q11, q12, q5, q6), ncol = 2, byrow = FALSE)
            partition <- matrix(data, nrow = 2, byrow = FALSE)
            C[, (2:nn)] <- QB %*% partition[, (2:nn) - 1] + QA %*% 
                partition[, (2:nn)] + QZ %*% partition[, (2:nn) +
                1]
            C[, 1] <- QB %*% partition[, nn + 1] + QA %*% partition[
                , 1] + QZ %*% partition[, 2]
            C[, nn + 1] <- QB %*% partition[, nn] + QA %*% 
                partition[, nn + 1] + QZ %*% partition[, 1]
        }
        else stop("Bad prefilter for specified multiwavelet filter")
    }
    else if(filter.type == "Donovan3") {
        if(prefilter.type == "Identity") {
            if(verbose == TRUE)
                cat("  O.K.\nPrefilter (Identity)...")
            for(l in 1:nphi) {
                C[l, 1:(ndata/nphi)] <- data[(0:((ndata/nphi) - 
                  1)) * nphi + l]
            }
        }
        else if(prefilter.type == "Linear") {
            if(verbose == TRUE)
                cat("  O.K.\nPrefilter (Linear)...")
            C[1, 1:(ndata/3)] <- data[3 * 0:((ndata/3) - 1) + 1] * 
                -0.76885512
                
 + data[3 * 0:((ndata/3) - 1) + 2]
            C[2, 1:(ndata/3)] <- data[3 * 0:((ndata/3) - 1) + 1] * 
                -0.56536682999999999
                
 + data[3 * 0:((ndata/3) - 1) + 2]
            C[3, 1:(ndata/3)] <- data[3 * 0:((ndata/3) - 1) + 1] * 
                0.097676540000000006 - data[3 * 0:((ndata/3) - 
                1) + 2] + data[3 * 1:(ndata/3)]
        }
        else if(prefilter.type == "Interp" || prefilter.type == 
            "default") {
            if(verbose == TRUE)
                cat("  O.K.\nPrefilter (Interpolation)...")
            w <- sqrt(5)
            lc <- length(data)/3
            C[1, 1:lc] <- data[3 * (0:(lc - 1)) + 1] * sqrt(7/11)
            C[3, 1] <- ((sqrt(3) * (data[2] - data[3]) + (C[1, lc] * (
                -1 + 8 * w))/3/sqrt(231) + (C[1, 1] * (1 + 8 * 
                w))/3/sqrt(231)) * 3 * sqrt(33) * (16 - 5 * w))/
                (-203 + 88 * w)
            C[3, 2:lc] <- ((sqrt(3) * (data[3 * (1:(lc - 1)) + 2] - 
                data[3 * (2:lc)]) + (C[1, 1:(lc - 1)] * (-1 + 8 *
                w))/3/sqrt(231) + (C[1, 2:lc] * (1 + 8 * w))/3/
                sqrt(231)) * 3 * sqrt(33) * (16 - 5 * w))/(-203 +
                88 * w)
            C[2, 1] <- ((sqrt(3) * data[2] + (C[1, lc] * (2 + 6 * w
                ))/3/sqrt(231) + (C[1, 1] * (3 + 2 * w))/3/sqrt(
                231) - (C[3, 1] * (103 - 24 * w))/3/sqrt(33)/(
                16 - 5 * w)) * sqrt(3))/2
            C[2, 2:lc] <- ((sqrt(3) * data[3 * (1:(lc - 1)) + 2] + (
                C[1, 1:(lc - 1)] * (2 + 6 * w))/3/sqrt(231) + (
                C[1, 2:lc] * (3 + 2 * w))/3/sqrt(231) - (C[3, 2:
                lc] * (103 - 24 * w))/3/sqrt(33)/(16 - 5 * w)) * 
                sqrt(3))/2
        }
        else stop("Bad prefilter for specified multiwavelet filter")
    }
    else stop("No prefilter for the multiwavelet filter")
    return(C)
}
"mwd"<-
function(data, prefilter.type = "default", filter.type = "Geronimo", bc = 
    "periodic", verbose = FALSE)
{
#
#applies the Discrete Multiple wavelet Transform to data
#copyrigt Tim Downie 1995-1996
#
    if(verbose == TRUE) cat("Multiple wavelet decomposition\n")
    if(verbose == TRUE)
        cat("Checking Arguements...")
    if(bc != "periodic")
        stop("\nOnly periodic boundary conditions allowed at the moment"
            )
    filter <- mfilter.select(type = filter.type)
    ndata <- length(data)   #   
#
# check ndata = filter$nphi * filter$ndecim ^ nlevels 
# 
#
    nlevels <- log(ndata/filter$nphi)/log(filter$ndecim)    #
#
#  repeated signal prefilter has one extra level
#
    if(prefilter.type == "Repeat")
        nlevels <- nlevels + 1
    if(nlevels != round(nlevels) || nlevels < 1)
        stop("\nbad number of data points for this filter\n")
    if(verbose == TRUE)
        cat("  O.K.\nBuilding first/last database ...")
    fl <- mfirst.last(LengthH = filter$NH, nlevels = nlevels, ndecim = 
        filter$ndecim, type = "wavelet", bc = bc)   #
    if(bc == "periodic")
        nbc <- 1
    else if(bc == "symmetric")
        nbc <- 2
    C <- mprefilter(data, prefilter.type, filter.type, nlevels, fl$nvecs.c, 
        filter$nphi, filter$npsi, filter$ndecim, verbose)
    if(verbose == TRUE)
        cat(" O.K.\nRunning decomposition algorithm...")
    gwd <- .C("multiwd",
        C = as.double(C),
        lengthc = as.integer(fl$nvecs.c * filter$nphi),
        D = as.double(rep(0, fl$nvecs.d * filter$npsi)),
        lengthd = as.integer(fl$nvecs.d * filter$npsi),
        nlevels = as.integer(nlevels),
        nphi = as.integer(filter$nphi),
        npsi = as.integer(filter$npsi),
        ndecim = as.integer(filter$ndecim),
        H = as.double(filter$H),
        G = as.double(filter$G),
        NH = as.integer(filter$NH),
        lowerc = as.integer(fl$first.last.c[, 1]),
        upperc = as.integer(fl$first.last.c[, 2]),
        offsetc = as.integer(fl$first.last.c[, 3]),
        lowerd = as.integer(fl$first.last.d[, 1]),
        upperd = as.integer(fl$first.last.d[, 2]),
        offsetd = as.integer(fl$first.last.d[, 3]),
        nbc = as.integer(nbc), PACKAGE = "wavethresh")  # 
# the C function returns the C and D coefficients as a vector
# convert into a matrix with nphi rows.
# 
    gwd$C <- matrix(gwd$C, nrow = filter$nphi)
    gwd$D <- matrix(gwd$D, nrow = filter$npsi)
    outlist <- list(C = gwd$C, D = gwd$D, nlevels = nlevels, ndata = ndata, 
        filter = filter, fl.dbase = fl, type = "wavelet", bc = bc, 
        prefilter = prefilter.type, date = date())
    class(outlist) <- "mwd"
    if(verbose == TRUE)
        cat(" O.K.\nReturning Multiple Wavelet Decomposition\n")
    return(outlist)
}
"mwr"<-
function(mwd, prefilter.type = mwd$prefilter, verbose = FALSE, start.level = 0, 
    returnC = FALSE)
{
#function to reconstruct the data from an object of class mwd
#a multiwavelet decomposition
#Tim Downie
#last updated May 96
    if(verbose == TRUE) cat("Multiple wavelet reconstruction\nArguement checking ..."
            )
    ctmp <- class(mwd)
    if(is.null(ctmp))
        stop("Input must have class mwd")
    else if(ctmp != "mwd")
        stop("Input must have class mwd")
    if(mwd$prefilter != prefilter.type)
        warning("The pre/postfilters are inconsistent\n")
    if(start.level < 0 || start.level >= nlevelsWT(mwd)) stop(
            "Start.level out of range\n")   #
# keep the value of the Cs at level 0 reset all the others
#
    if(verbose == TRUE)
        cat(" O.K.\nInitialising variables ...")
    C <- matrix(rep(0, length(mwd$C)), nrow = mwd$filter$nphi)
    c0low <- mwd$fl.dbase$first.last.c[start.level + 1, 3] + 1
    c0high <- c0low + mwd$fl.dbase$first.last.c[start.level + 1, 2] - mwd$
        fl.dbase$first.last.c[start.level + 1, 1]
    for(l in 1:mwd$filter$nphi)
        C[l, c0low:c0high] <- mwd$C[l, c0low:c0high]
    if(mwd$bc == "periodic")
        nbc <- 1
    else if(mwd$bc == "symmetric")
        nbc <- 2
    else stop("bad boundary conditions")
    if(verbose == TRUE)
        cat(" O.K.\nRunning Reconstruction algorithm...")
    reconstr <- .C("multiwr",
        C = as.double(C),
        lengthc = as.integer(mwd$fl.dbase$ntotal),
        D = as.double(mwd$D),
        lengthd = as.integer(mwd$fl.dbase$ntotal.d),
        nlevels = as.integer(nlevelsWT(mwd)),
        nphi = as.integer(mwd$filter$nphi),
        npsi = as.integer(mwd$filter$npsi),
        ndecim = as.integer(mwd$filter$ndecim),
        H = as.double(mwd$filter$H),
        G = as.double(mwd$filter$G),
        NH = as.integer(mwd$filter$NH),
        lowerc = as.integer(mwd$fl.dbase$first.last.c[, 1]),
        upperc = as.integer(mwd$fl.dbase$first.last.c[, 2]),
        offsetc = as.integer(mwd$fl.dbase$first.last.c[, 3]),
        lowerd = as.integer(mwd$fl.dbase$first.last.d[, 1]),
        upperd = as.integer(mwd$fl.dbase$first.last.d[, 2]),
        offsetd = as.integer(mwd$fl.dbase$first.last.d[, 3]),
        nbc = as.integer(nbc),
        startlevel = as.integer(start.level), PACKAGE = "wavethresh")
    ndata <- mwd$filter$ndecim^nlevelsWT(mwd)* mwd$filter$nphi
    reconstr$C <- matrix(reconstr$C, nrow = mwd$filter$nphi)
    if(returnC == TRUE) {
        if(verbose == TRUE)
            cat(" O.K.\nReturning starting coefficients\n")
        return(reconstr$C[, (1:(ndata/mwd$filter$nphi))])
    }
    if(verbose == TRUE)
        cat(" O.K.\nApply post filter...")
    ndata <- mwd$filter$ndecim^nlevelsWT(mwd)* mwd$filter$nphi
    data <- mpostfilter(reconstr$C, prefilter.type, mwd$filter$type, mwd$
        filter$nphi, mwd$filter$npsi, mwd$filter$ndecim, nlevelsWT(mwd), 
        verbose)
    if(verbose == TRUE)
        cat(" O.K.\nReturning data\n")
    return(data)
}
"newsure"<-
function(s, x)
{
    x <- abs(x)
    d <- length(x)
    sl <- sort.list(x)
    y <- x[sl]
    sigma <- s[sl]
    cy <- cumsum(y^2)
    cy <- c(0, cy[1:(length(cy) - 1)])
    csigma <- cumsum(sigma^2)
    csigma <- c(0, csigma[1:(length(csigma) - 1)])
    ans <- d - 2 * csigma + cy + d:1 * y^2
    m <- min(ans)
    index <- (1:length(ans))[m == ans]
    return(y[index])
}
"nlevelsWT"<-
function(...)
UseMethod("nlevelsWT")

#"nlevels.default"<-
#function(object, ...)
#{
#    if(is.null(object$nlevels)) {
#        n <- length(object)
#        return(IsPowerOfTwo(n))
#    }
#    else return(object$nlevels)
#}

#MAN: changed function below to cope with $nlevels deprecation (R-2.6.0 onwards).

"nlevelsWT.default"<-
function(object, ...)
{
if (is.list(object)){
    if(!is.null(object$nlevels)){       # "normal" object */
        return(object$nlevels)
    }
    else{
        if(class(object)=="uncompressed"){      # 2 special cases 
            return(IsPowerOfTwo(object$v))
        }
        else if(class(object)=="griddata"){
            return(IsPowerOfTwo(object$gridy))
        }
        else{                                       # what to do?  e.g. tpwd,wpstDO,compressed classes. 
            print("I don't know what to do with this object!\n")
            stop("unknown nlevels")
        }

    }
}    
else{                                           #data should be atomic (numeric)...
        return(IsPowerOfTwo(length(object)))
}

}


"nullevels"<-
function(...)
UseMethod("nullevels")
"nullevels.imwd"<-
function(imwd, levelstonull, ...)
{
    nlevels <- nlevelsWT(imwd)
    if(max(levelstonull) > nlevels - 1)
        stop(paste("Illegal level to null, maximum is ", nlevels - 1))
    if(min(levelstonull) < 0)
        stop(paste("Illegal level to null, minimum is ", nlevels - 1))
    for(lev in levelstonull) {
        n1 <- lt.to.name(lev, type = "CD")
        n2 <- lt.to.name(lev, type = "DC")
        n3 <- lt.to.name(lev, type = "DD")
        imwd[[n1]] <- rep(0, length(imwd[[n1]]))
        imwd[[n2]] <- rep(0, length(imwd[[n2]]))
        imwd[[n3]] <- rep(0, length(imwd[[n3]]))
    }
    imwd
}
"nullevels.wd"<-
function(wd, levelstonull, ...)
{
    nlevels <- nlevelsWT(wd)
    if(max(levelstonull) > nlevels - 1)
        stop(paste("Illegal level to null, maximum is ", nlevels - 1))
    if(min(levelstonull) < 0)
        stop(paste("Illegal level to null, minimum is ", nlevels - 1))
    for(lev in levelstonull) {
        d <- accessD(wd, level = lev)
        d <- rep(0, length(d))
        wd <- putD(wd, level = lev, v = d)
    }
    wd
}
"nullevels.wst"<-
function(wst, levelstonull, ...)
{
    nullevels.wd(wst, levelstonull = levelstonull)
}
"numtonv"<-
function(number, nlevels)
{
    if(nlevels < 1)
        stop("nlevels cannot be less than 1")
    if(number < 0)
        stop("Number cannot be less than 0")
    else if(number > 2^nlevels - 1)
        stop(paste("Number cannot be more than", 2^nlevels - 1))
    node.vector <- vector("list", nlevels)
    matchcodes <- c("L", "R")
    mask <- 2^(nlevels - 1)
    cmc <- NULL
    for(i in (nlevels - 1):0) {
        index <- floor(number/mask)
        if(index == 1)
            number <- number - mask
        mask <- mask/2
        cmc <- c(cmc, index)
    }
    for(i in (nlevels - 1):0) {
        index <- cmc[i + 1]
        nul <- 2^(nlevels - i - 1)
        upperl <- rep(0, nul)
        upperctrl <- rep(matchcodes[index + 1], nul)
        node.vector[[i + 1]] <- list(upperctrl = upperctrl, upperl = 
            upperl)
    }
    node.vector <- list(node.list = node.vector, nlevels = nlevels)
    class(node.vector) <- "nv"
    node.vector
}

"plot.imwd"<-
function(x, scaling = "by.level", co.type = "abs", package = "R", 
    plot.type = "mallat", arrangement = c(3, 3), transform = FALSE, tfunction
     = sqrt, ...)
{
#
#
#       Check class of imwd
#
    if(package != "R" && package != "S") stop("Unknown package")
    ctmp <- class(x)
    if(is.null(ctmp))
        stop("imwd has no class")
    else if(ctmp != "imwd")
        stop("imwd is not of class imwd")
    if(x$type == "station" && plot.type == "mallat")
        stop("Cannot do Mallat type plot on nondecimated wavelet object")
    Csize <- 2^(nlevelsWT(x))
    m <- matrix(0, nrow = Csize, ncol = Csize)
    first.last.d <- x$fl.dbase$first.last.d
    first.last.c <- x$fl.dbase$first.last.c
    if(plot.type == "mallat") {
        for(level in (nlevelsWT(x)):1) {
            ndata <- 2^(level - 1)
            firstD <- first.last.d[level, 1]
            lastD <- first.last.d[level, 2]
            LengthD <- lastD - firstD + 1
            sel <- seq(from = (1 - firstD), length = ndata) #
#
# Extract CD for this level
#
            nm <- lt.to.name(level - 1, "CD")
            msub1 <- matrix(x[[nm]], nrow = LengthD, ncol = 
                LengthD)    #
#
# Extract DC for this level
#
            nm <- lt.to.name(level - 1, "DC")
            msub2 <- matrix(x[[nm]], nrow = LengthD, ncol = 
                LengthD)    #
#
# Extract DD for this level
#
            nm <- lt.to.name(level - 1, "DD")
            msub3 <- matrix(x[[nm]], nrow = LengthD, ncol = 
                LengthD)    #
#
#
#   Work out if we want to display the absolute values or the actual
#   values
#
            if(co.type == "abs") {
                msub1 <- abs(msub1)
                msub2 <- abs(msub2)
                msub3 <- abs(msub3)
            }
            else if(co.type == "mabs") {
                msub1 <-  - abs(msub1)
                msub2 <-  - abs(msub2)
                msub3 <-  - abs(msub3)
            }
            else if(co.type != "none")
                stop("Unknown co.type")
            if(transform == TRUE) {
                msub1 <- tfunction(msub1)
                msub2 <- tfunction(msub2)
                msub3 <- tfunction(msub3)
            }
            if(scaling == "by.level") {
                if(ndata == 1) {
                  r.m1 <- range(c(as.vector(msub1), as.vector(
                    msub2), as.vector(msub3)))
                  r.m2 <- r.m1
                  r.m3 <- r.m1
                }
                else {
                  r.m1 <- range(msub1)
                  r.m2 <- range(msub2)
                  r.m3 <- range(msub3)
                }
                if(r.m1[2] - r.m1[1] == 0) {
                  msub1[,  ] <- 0
                }
                else {
                  mu1 <- 249/(r.m1[2] - r.m1[1])
                  msub1 <- mu1 * (msub1 - r.m1[1])
                }
                if(r.m2[2] - r.m2[1] == 0) {
                  msub2[,  ] <- 0
                }
                else {
                  mu2 <- 249/(r.m2[2] - r.m2[1])
                  msub2 <- mu2 * (msub2 - r.m2[1])
                }
                if(r.m3[2] - r.m3[1] == 0) {
                  msub3[,  ] <- 0
                }
                else {
                  mu3 <- 249/(r.m3[2] - r.m3[1])
                  msub3 <- mu3 * (msub3 - r.m3[1])
                }
            }
            else {
                range.msub <- range(c(msub1, msub2, msub3))
                multiplier <- 255/(range.msub[2] - range.msub[1
                  ])
                msub1 <- multiplier * (msub1 - range.msub[1])
                msub2 <- multiplier * (msub2 - range.msub[1])
                msub3 <- multiplier * (msub3 - range.msub[1])   #
            }
            m[(ndata + 1):(2 * ndata), 1:ndata] <- msub1[sel, sel]
            m[1:ndata, (ndata + 1):(2 * ndata)] <- msub2[sel, sel]
            m[(ndata + 1):(2 * ndata), (ndata + 1):(2 * ndata)] <- 
                msub3[sel, sel]
        }
        if(package == "R") {
            image(m, xaxt = "n", yaxt = "n",...)
            axis(1, at = c(0, 2^((nlevelsWT(x)- 3):(nlevelsWT(x)))
                ))
            axis(2, at = c(0, 2^((nlevelsWT(x)- 3):(nlevelsWT(x)))
                ))
        }
        else return(m)
    }
    else if(plot.type == "cols") {
        oldpar <- par(mfrow = arrangement, pty = "s")
        for(level in (nlevelsWT(x):1)) {
            ndata <- 2^(level - 1)
            firstD <- first.last.d[level, 1]
            lastD <- first.last.d[level, 2]
            LengthD <- lastD - firstD + 1
            sel <- seq(from = (1 - firstD), length = ndata) #
#
# Extract CD for this level
#
            nm <- lt.to.name(level - 1, "CD")
            msub1 <- matrix(x[[nm]], nrow = LengthD, ncol = 
                LengthD)    #
#
# Extract DC for this level
#
            nm <- lt.to.name(level - 1, "DC")
            msub2 <- matrix(x[[nm]], nrow = LengthD, ncol = 
                LengthD)    #
#
# Extract DD for this level
#
            nm <- lt.to.name(level - 1, "DD")
            msub3 <- matrix(x[[nm]], nrow = LengthD, ncol = 
                LengthD)    #
#
#
#   Work out if we want to display the absolute values or the actual
#   values
#
            if(co.type == "abs") {
                msub1 <- abs(msub1)
                msub2 <- abs(msub2)
                msub3 <- abs(msub3)
            }
            else if(co.type == "mabs") {
                msub1 <-  - abs(msub1)
                msub2 <-  - abs(msub2)
                msub3 <-  - abs(msub3)
            }
            else if(co.type != "none")
                stop("Unknown co.type")
            if(transform == TRUE) {
                msub1 <- tfunction(msub1)
                msub2 <- tfunction(msub2)
                msub3 <- tfunction(msub3)
            }
            if(package == "R") {
                xlabstr <- paste("Level", level - 1, 
                  "(horizonatal)")
                image(msub1, xlab = xlabstr)
                xlabstr <- paste("Level", level - 1, 
                  "(vertical)")
                image(msub2, xlab = xlabstr)
                xlabstr <- paste("Level", level - 1, 
                  "(diagonal)")
                image(msub3, xlab = xlabstr,...)
            }
            else {
                warning("Not using R")
            }
        }
        par(oldpar)
    }
    else stop("Unknown plot.type")
}
"plot.imwdc"<-
function(x, verbose = FALSE, ...)
{
    imwd <- uncompress(x, verbose = verbose)
    return(plot(imwd, ...))
}
plot.irregwd <-
function (x, xlabels, first.level = 1, main = "Wavelet Decomposition Coefficients", 
    scaling = "by.level", rhlab = FALSE, sub, ...) 
{
    ctmp <- class(x)
    if (is.null(ctmp)) 
        stop("irregwd has no class")
    else if (ctmp != "irregwd") 
        stop("irregwd is not of class irregwd")
    iwd <- x
    wd <- x
    class(wd) <- "wd"
    levels <- nlevelsWT(wd)
    nlevels <- levels - first.level
    n <- 2^(levels - 1)
    if (missing(sub)) 
        sub <- wd$filter$name
    plot(c(0, 0, n, n), c(0, nlevels + 1, nlevels + 1, 0), type = "n", 
        xlab = "Translate", ylab = "Resolution Level", main = main, 
        yaxt = "n", xaxt = "n", sub = sub, ...)
    axis(2, at = 1:(nlevels), labels = ((levels - 1):first.level))
    if (missing(xlabels)) {
        axx <- c(0, 2^(nlevels - 2), 2^(nlevels - 1), 2^(nlevels - 
            1) + 2^(nlevels - 2), 2^nlevels)
        axis(1, at = axx)
    }
    else {
        axx <- pretty(1:n, n = 3)
        if (axx[length(axx)] > n) 
            axx[length(axx)] <- n
        axx[axx == 0] <- 1
        axl <- signif(xlabels[axx], digits = 3)
        axis(1, at = axx, labels = axl)
    }
    x <- 1:n
    height <- 1
    first.last.d <- wd$fl.dbase$first.last.d
    axr <- NULL
    if (scaling == "global") {
        my <- 0
        for (i in ((levels - 1):first.level)) {
            y <- accessc(iwd, i)
            my <- max(c(my, abs(y)))
        }
    }
    for (i in ((levels - 1):first.level)) {
        n <- 2^i
        y <- accessc(iwd, i)
        xplot <- x
        ly <- length(y)
        if (scaling == "by.level") 
            my <- max(abs(y))
        y <- (0.5 * y)/my
        axr <- c(axr, my)
        segments(xplot, height, xplot, height + y)
        if (i != first.level) {
            x1 <- x[seq(1, n - 1, 2)]
            x2 <- x[seq(2, n, 2)]
            x <- (x1 + x2)/2
            height <- height + 1
        }
    }
    if (rhlab == TRUE) 
        axis(4, at = 1:length(axr), labels = signif(axr, 3))
    axr
}

"plot.mwd"<-
function(x, first.level = 1, main = "Wavelet Decomposition Coefficients", 
    scaling = "compensated", rhlab = FALSE, sub = x$filter$name, NotPlotVal
     = 0.050000000000000003, xlab = "Translate", ylab = "Resolution level", 
    return.scale = TRUE, colour = (2:(npsi + 1)), ...)
{
#plot.mwd
#plot a multiwavelet decompostion
#
#Tim Downie  1995-1996
#
#
#       Check class of mwd
#
    ctmp <- class(x)
    if(is.null(ctmp))
        stop("mwd has no class")
    else if(ctmp == "wd")
        stop("object is of class wd use plot.wd or plot")
    else if(ctmp != "mwd")
        stop("object is not of class mwd")
    nlevels <- nlevelsWT(x)- first.level
    mx <- x$ndata
    xlabs <- seq(0, mx/2, length = 5)
    plot(c(0, 0, mx, mx), c(0, nlevels + 1, nlevels + 1, 0), type = "n", 
        xlab = xlab, ylab = ylab, main = main, yaxt = "n", xaxt = "n", 
        sub=sub, ...)
    axis(1, at = seq(0, mx, length = 5), labels = xlabs)
    axis(2, at = 1:(nlevels), labels = (nlevelsWT(x)- 1):first.level)
    delta <- 1
    npsi <- x$filter$npsi
    ndecim <- x$filter$ndecim
    height <- 1
    first.last.d <- x$fl.dbase$first.last.d
    axr <- NULL
    if(scaling == "global") {
        my <- 0
        for(i in ((nlevelsWT(x)- 1):first.level)) {
            y <- c(accessD(x, i))
            my <- max(c(my, abs(y)))
        }
    }
    if(scaling == "compensated") {
        my <- 0
        for(i in ((nlevelsWT(x)- 1):first.level)) {
            y <- c(accessD(x, i)) * x$filter$ndecim^(i/2)
            my <- max(c(my, abs(y)))
        }
    }
    for(i in ((nlevelsWT(x)- 1):first.level)) {
        y <- c(accessD(x, i))
        ly <- length(y)
        n <- ly/npsi
        if(scaling == "by.level")
            my <- max(abs(y))
        if(scaling == "compensated")
            y <- y * ndecim^(i/2)
        if(my == 0)
            y <- rep(0, ly)
        else y <- (0.5 * y)/my
        axr <- c(axr, my)
        xplot <- rep(((1:n) * mx)/(n + 1), rep(npsi, ly/npsi)) + (0:(
            npsi - 1)) * delta
        segments(xplot, height, xplot, height + y, col = colour)
        height <- height + 1
    }
    if(rhlab == TRUE)
        axis(4, at = 1:length(axr), labels = signif(axr, 3))
    if(return.scale == TRUE)
        return(axr)
    else return(NULL)
}
"plot.nvwp"<-
function(x, ...)
{
    plotpkt(nlevelsWT(x))
    pktlist <- print.nvwp(x, printing = FALSE)
    for(i in 1:length(pktlist$level))
        addpkt(pktlist$level[i], pktlist$pkt[i], 1, col = 1)
}
"plot.wd"<-
function(x, xlabvals, xlabchars, ylabchars, first.level = 0, main = 
    "Wavelet Decomposition Coefficients", scaling = "global", rhlab = FALSE, 
    sub, NotPlotVal = 0.0050000000000000001, xlab = "Translate", ylab = 
    "Resolution Level", aspect = "Identity", ...)
{
    if(IsEarly(x)) {
        ConvertMessage()
        stop()
    }
    if(is.complex(x$D) && aspect == "Identity") aspect <- "Mod"    #
#       Check class of wd
#
    ctmp <- class(x)
    if(is.null(ctmp))
        stop("wd has no class")
    else if(ctmp != "wd")
        stop("wd is not of class wd")
    levels <- nlevelsWT(x)
    if(x$bc == "interval") {
        if(first.level < x$current.scale)
            warning(paste("plot.wd plotted from level", x$
                current.scale, 
                " because \"wavelets on the interval\" transform was only computed to this level\n"
                ))
        first.level <- x$current.scale
    }
    nlevels <- levels - first.level
    type <- x$type
    if(IsEarly(x)) {
        ConvertMessage()
        stop()
    }
    if(type == "wavelet")
        n <- 2^(levels - 1)
    else if(type == "station")
        n <- 2^levels
    else stop("Unknown type for wavelet object")
    if(missing(sub))
        sub <- paste(switch(type,
            wavelet = "Standard transform",
            station = "Nondecimated transform"), x$filter$name)
    if(aspect != "Identity")
        sub <- paste(sub, "(", aspect, ")")
    plot(c(0, 0, n, n), c(0, nlevels + 1, nlevels + 1, 0), type = "n", xlab
         = xlab, ylab = ylab, main = main, yaxt = "n", xaxt = "n", sub
         = sub, ...)
    yll <- (levels - 1):first.level
    if(missing(ylabchars))
        axis(2, at = 1:(nlevels), labels = yll)
    else if(length(ylabchars) != nlevels)
        stop(paste("Should have ", nlevels, " entries in ylabchars"))
    else axis(2, at = 1:(nlevels), labels = ylabchars)
    if(missing(xlabchars)) {
        if(missing(xlabvals)) {
            if(type == "wavelet")
                axx <- c(0, 2^(levels - 3), 2^(levels - 2), 2^(
                  levels - 2) + 2^(levels - 3), 2^(levels - 1))
            else axx <- c(0, 2^(levels - 2), 2^(levels - 1), 2^(
                  levels - 1) + 2^(levels - 2), 2^levels)
            if(is.null(tsp(x)))
                axis(1, at = axx)
            else {
                v <- seq(from = tsp(x)["start"], by = tsp(
                  x)["deltat"], length = n)
                if(type == "wavelet")
                  atl <- 2 * v
                else atl <- v
                atl <- pretty(atl, n = 4)
                ats <- (n * atl)/(max(atl) - min(atl))
                axis(1, at = ats, labels = atl)
            }
        }
        else {
            lx <- pretty(xlabvals, n = 4)
            cat("lx is ", lx, "\n")
            if(lx[1] < min(xlabvals))
                lx[1] <- min(xlabvals)
            if(lx[length(lx)] > max(xlabvals))
                lx[length(lx)] <- max(xlabvals)
            cat("lx is ", lx, "\n")
            xix <- NULL
            for(i in 1:length(lx)) {
                u <- (xlabvals - lx[i])^2
                xix <- c(xix, (1:length(u))[u == min(u)])
            }
            axx <- xix
            if(type == "wavelet")
                axx <- xix/2
            axl <- signif(lx, digits = 2)
            axis(1, at = axx, labels = axl)
        }
    }
    else axis(1, at = xlabvals, labels = xlabchars)
    myxx <- 1:n
    height <- 1
    first.last.d <- x$fl.dbase$first.last.d
    axr <- NULL
    if(scaling == "global") {
        my <- 0
        for(i in ((levels - 1):first.level)) {
            y <- accessD(x, i, aspect = aspect)
            my <- max(c(my, abs(y)))
        }
    }
    if(scaling == "compensated") {
        my <- 0
        for(i in ((levels - 1):first.level)) {
            y <- accessD(x, i, aspect = aspect) * 2^(i/2)
            my <- max(c(my, abs(y)))
        }
    }
    if(scaling == "super") {
        my <- 0
        for(i in ((levels - 1):first.level)) {
            y <- accessD(x, i, aspect = aspect) * 2^i
            my <- max(c(my, abs(y)))
        }
    }
    shift <- 1
    for(i in ((levels - 1):first.level)) {
        y <- accessD(x, i, aspect = aspect)
        if(type == "wavelet")
            n <- 2^i
        else {
            y <- y[c((n - shift + 1):n, 1:(n - shift))]
            shift <- shift * 2
        }
        xplot <- myxx
        ly <- length(y)
        if(scaling == "by.level")
            my <- max(abs(y))
        if(scaling == "compensated")
            y <- y * 2^(i/2)
        if(scaling == "super")
            y <- y * 2^i
        if(my == 0) {
            y <- rep(0, length(y))
        }
        else y <- (0.5 * y)/my
        axr <- c(axr, my)
        if(max(abs(y)) > NotPlotVal)
            segments(xplot, height, xplot, height + y)
        if(i != first.level) {
            if(type == "wavelet") {
                x1 <- myxx[seq(1, n - 1, 2)]
                x2 <- myxx[seq(2, n, 2)]
                myxx <- (x1 + x2)/2
            }
            height <- height + 1
        }
    }
    if(rhlab == TRUE)
        axis(4, at = 1:length(axr), labels = signif(axr, digits=3))
    axr
}
"plot.wp"<-
function(x, nvwp = NULL, main = "Wavelet Packet Decomposition", sub, 
    first.level = 5, scaling = "compensated", dotted.turn.on = 5, 
    color.force = FALSE, WaveletColor = 2, NodeVecColor = 3, fast = FALSE, 
    SmoothedLines = TRUE, ...)
{
#
# Check class of wp
#
    ctmp <- class(x)
    if(is.null(ctmp))
        stop("wp has no class")
    else if(ctmp != "wp")
        stop("wp is not of class wp")
    levels <- nlevelsWT(x)
    dotted.turn.on <- levels - dotted.turn.on
    N <- 2^levels   # The number of original data points
#
#
#   Check validity of command line args
#
    if(first.level < 0 || first.level > levels)
        stop("first.level must between zero and the number of levels")  
    #
    if(dotted.turn.on < 0 || dotted.turn.on > levels) stop(
            "dotted.turn.on must between zero and number of levels"
            )   #
#   Do subtitling
#
    if(missing(sub)) sub <- paste("Filter: ", x$filter$name)   #
#
#   Set plotting region and do axes of plot
#
    oldpar <- par(mfrow = c(1, 1))
    if(!is.null(nvwp))
        sub <- paste(sub, "(selected packets in color 3)")
    plot(c(0, N + 1), c(-1, levels - first.level + 1), type = "n", main = 
        main, xlab = "Packet Number", ylab = "Resolution Level", yaxt
         = "n", sub = sub, ...)
    axis(2, at = 0:(levels - first.level), labels = levels:first.level) #
#
#   Check out how to do things in a different colour if we can
#
    if(color.force == FALSE) {
        if(CanUseMoreThanOneColor() == FALSE) {
            if(WaveletColor > 1) {
                warning(
                  "Can't (or can't find out how) display wavelets in color"
                  )
                WaveletColor <- 1
            }
            if(NodeVecColor > 1) {
                warning(
                  "Can't (or can't find out how) display node vector packets in color"
                  )
                NodeVecColor <- 1
            }
        }
    }
    origdata <- getpacket(x, lev = levels, index = 0)  #
#
#   Scaling for the original data is always the same
#
    sf <- max(abs(origdata))
    if(sf == 0) {
        stop("Original data is the zero function\n")
    }
    scale.origdata <- (0.5 * origdata)/sf
    lines(1:N, scale.origdata)
    if(first.level == levels) return()  #
#
#   Draw the vertical seperators if necessary
#
    for(i in 1:(levels - first.level)) {
        N <- N/2
        if(i > dotted.turn.on)
            break
        else for(j in 1:(2^i - 1)) {
                segments(N * (j - 1) + N + 0.5, i - 0.5, N * (j -
                  1) + N + 0.5, i + 0.5, lty = 2)
            }
    }
#
#
#   Get all the coefficients    
#
    CoefMatrix <- x$wp #
#
#   Remove the original data cos we've already plotted that
#
    CoefMatrix <- CoefMatrix[ - (levels + 1),  ]    #
#   Compute Global Scale Factor if necessary
#
    Sf <- 0
    if(scaling == "global")
        Sf <- max(abs(CoefMatrix), na.rm = TRUE)
    else if(scaling == "compensated") {
        for(i in 1:(levels - first.level)) {
            Coefs <- CoefMatrix[levels - i + 1,  ] * 2^((levels - i
                )/2)
            Sf <- max(c(Sf, abs(Coefs)), na.rm = TRUE)
        }
    }
    if(scaling == "global")
        sf <- Sf
    if(is.null(nvwp)) {
#
#   If there is no associated node vector then plot the wavelet packet
#   table using the matrix of coefficients. This is faster than the
#   packet by packet method that is used when we have a node vector
#   (but probably not much)
#
#
        for(i in 1:(levels - first.level)) {
            PKLength <- 2^(levels - i)
            Coefs <- CoefMatrix[levels - i + 1,  ]
            if(scaling == "by.level")
                sf <- max(abs(Coefs), na.rm = TRUE)
            else if(scaling == "compensated")
                sf <- Sf/2^((levels - i)/2)
            if(is.na(sf) || sf == 0)
                Coefs <- rep(0, length(Coefs))
            else Coefs <- (0.5 * Coefs)/sf
            pkl <- 1:PKLength
            if(SmoothedLines == TRUE)
                lines(pkl, i + Coefs[pkl])
            else segments(pkl, i, pkl, i + Coefs[pkl])
            pkl <- PKLength + pkl
            segments(pkl, i, pkl, i + Coefs[pkl], col=WaveletColor)
            pkl <- (2 * PKLength + 1):length(Coefs)
            segments(pkl, i, pkl, i + Coefs[pkl])
        }
    }
    else {
        pklist <- print.nvwp(nvwp, printing = FALSE)
        for(i in 1:(levels - first.level)) {
#
#           Scaling issues
#
            Coefs <- CoefMatrix[levels - i + 1,  ]
            if(scaling == "by.level")
                sf <- max(abs(Coefs), na.rm = TRUE)
            else if(scaling == "compensated")
                sf <- Sf/2^((levels - i)/2)
            if(is.na(sf) || sf == 0)
                Coefs <- rep(0, length(Coefs))
            else Coefs <- (0.5 * Coefs)/sf
            CoefMatrix[levels - i + 1,  ] <- Coefs
            x$wp <- CoefMatrix
            the.lev <- levels - i
            PKLength <- 2^the.lev
            npkts <- 2^i
            pkl <- 1:PKLength
            for(j in 1:npkts) {
                pkt <- getpacket(x, level = the.lev, index = j -
                  1)
                lcol <- 1
                if(any(pklist$level == the.lev)) {
                  lpklist <- pklist$pkt[pklist$level == the.lev
                    ]
                  if(any(lpklist == (j - 1)))
                    lcol <- NodeVecColor
                  else if(j == 2)
                    lcol <- WaveletColor
                }
                else if(j == 2)
                  lcol <- WaveletColor
                if(j == 1) {
                  if(SmoothedLines == TRUE)
                    lines(pkl, i + pkt, col=lcol)
                  else segments(pkl, i, pkl, i + pkt, col=lcol)
                }
                else segments(pkl, i, pkl, i + pkt, col=lcol)
                pkl <- pkl + PKLength
            }
        }
    }
    invisible()
}
"plot.wst"<-
function(x, main = "Nondecimated Wavelet (Packet) Decomposition", sub, 
    first.level = 5, scaling = "compensated", dotted.turn.on = 5, aspect = 
    "Identity", ...)
{
#
# Check class of wst
#
    ctmp <- class(x)
    if(is.null(ctmp))
        stop("wst has no class")
    else if(ctmp != "wst")
        stop("wst is not of class wst")
    levels <- nlevelsWT(x)
    dotted.turn.on <- levels - dotted.turn.on
    if(is.complex(x$wp) && aspect == "Identity")
        aspect <- "Mod"
    N <- 2^levels   # The number of original data points
#
#
#   Check validity of command line args
#
    if(first.level < 0 || first.level > levels)
        stop("first.level must between zero and the number of levels")  
    #
    if(dotted.turn.on < 0 || dotted.turn.on > levels) stop(
            "dotted.turn.on must between zero and number of levels"
            )   #
#   Do subtitling
#
    if(missing(sub)) sub <- paste("Filter: ", x$filter$name)  #
#
#   Set plotting region and do axes of plot
#
    if(aspect != "Identity")
        sub <- paste(sub, "(", aspect, ")")
    plot(c(0, N + 1), c(-1, levels - first.level + 1), type = "n", main = 
        main, xlab = "Packet Number", ylab = "Resolution Level", yaxt
         = "n", sub = sub, ...)
    axis(2, at = 0:(levels - first.level), labels = levels:first.level) #
    origdata <- getpacket(x, lev = levels, index = 0, aspect = aspect)    #
#
#   Scaling for the original data is always the same
#
    sf <- max(abs(origdata))
    if(sf == 0) {
        scale.origdata <- rep(0, length(origdata))
    }
    else scale.origdata <- (0.5 * origdata)/sf
    lines(1:N, scale.origdata)
    if(first.level == levels) return()  #
#
#   Draw the vertical seperators if necessary
#
    for(i in 1:(levels - first.level)) {
        N <- N/2
        if(i > dotted.turn.on)
            break
        else for(j in 1:(2^i - 1)) {
                segments(N * (j - 1) + N + 0.5, i - 0.5, N * (j -
                  1) + N + 0.5, i + 0.5, lty = 2)
            }
    }
#
#
#   Get all the coefficients    
#
    if(aspect == "Identity")
        CoefMatrix <- x$wp
    else {
        fn <- get(aspect)
        CoefMatrix <- fn(x$wp)
    }
    CoefMatrix <- CoefMatrix[ - (levels + 1),  ]    #
#   Compute Global Scale Factor if necessary
#
    Sf <- 0
    if(scaling == "global")
        Sf <- max(abs(CoefMatrix), na.rm = TRUE)
    else if(scaling == "compensated") {
        for(i in 1:(levels - first.level)) {
            Coefs <- CoefMatrix[levels - i + 1,  ] * 2^((levels - i
                )/2)
            Sf <- max(c(Sf, abs(Coefs)), na.rm = TRUE)
        }
    }
    if(scaling == "global")
        sf <- Sf
    for(i in 1:(levels - first.level)) {
        PKLength <- 2^(levels - i)
        Coefs <- CoefMatrix[levels - i + 1,  ]
        if(scaling == "by.level")
            sf <- max(abs(Coefs), na.rm = TRUE)
        else if(scaling == "compensated")
            sf <- Sf/2^((levels - i)/2)
        if(is.na(sf) || sf == 0)
            Coefs <- rep(0, length(Coefs))
        else Coefs <- (0.5 * Coefs)/sf
        pkl <- 1:PKLength
        segments(pkl, i, pkl, i + Coefs[pkl])
        pkl <- PKLength + pkl
        segments(pkl, i, pkl, i + Coefs[pkl])
        pkl <- (2 * PKLength + 1):length(Coefs)
        segments(pkl, i, pkl, i + Coefs[pkl])
    }
}
"plot.wst2D"<-
function(x, plot.type = "level", main = "", ...)
{
    nlev <- nlevelsWT(x)
    sz <- dim(x$wst2D)[2]
    if(plot.type == "level") {
        for(i in 0:(nlev - 1)) {
            image(matrix(x$wst2D[i + 1,  ,  ], nrow = sz))
            st <- paste("Level", i)
            title(main = main, sub = st)
        }
    }
}
"plotpkt"<-
function(J)
{
    x <- c(0, 2^(J - 1))
    y <- c(0, J)
    plot(x, y, type = "n", xlab = "Packet indices", ylab = "Level", xaxt = 
        "n")
    axis(1, at = seq(from = 0, to = 2^(J - 1), by = 0.5), labels = 0:2^J)
}
"print.BP"<-
function(x, ...)
{
    cat("BP class object. Contains \"best basis\" information\n")
    cat("Components of object:")
    print(names(x))
    cat("Number of levels ", nlevelsWT(x), "\n")
    cat("List of \"best\" packets\n")
    m <- cbind(x$level, x$pkt, x$basiscoef)
    dimnames(m) <- list(NULL, c("Level id", "Packet id", "Basis coef"))
    print(m)
}
"print.imwd"<-
function(x, ...)
{
    cat("Class 'imwd' : Discrete Image Wavelet Transform Object:\n")
    cat("       ~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("$ wNLx are LONG coefficient vectors !\n")
    cat("\nsummary(.):\n----------\n")
    summary.imwd(x)
}
"print.imwdc"<-
function(x, ...)
{
    cat("Class 'imwdc' : Compressed Discrete Image Wavelet Transform Object:\n"
        )
    cat("       ~~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("$ wNLx are LONG coefficient vectors !\n")
    cat("\nsummary(.):\n----------\n")
    summary.imwdc(x)
}
"print.mwd"<-
function(x, ...)
{
    ctmp <- class(x)
    if(is.null(ctmp))
        stop("Input must have class mwd")
    else if(ctmp != "mwd")
        stop("Input must have class mwd")
    cat("Class 'mwd' : Discrete Multiple Wavelet Transform Object:\n")
    cat("       ~~~  : List with", length(x), "components with names\n")
    cat("              ", names(x), "\n\n")
    cat("$ C and $ D are LONG coefficient vectors !\n")
    cat("\nCreated on :", x$date, "\n")
    cat("Type of decomposition: ", x$type, "\n")
    cat("\nsummary:\n----------\n")
    summary.mwd(x)
}
"print.nv"<-
function(x, printing = TRUE, verbose = FALSE, ...)
{
    if(verbose == TRUE & printing == TRUE) {
        cat("Printing node vector as a list\n")
        cat("------------------------------\n")
        print(as.list(x))
        cat("Printing node vector as format\n")
        cat("------------------------------\n")
    }
    node.vector <- x$node.list
    acsel <- 0
    acsellist <- NULL
    cntr <- 0
    power <- 1
    rvector <- 0
    for(i in (nlevelsWT(x)- 1):0) {
        nl <- node.vector[[i + 1]]
        action <- nl$upperctrl[acsel + 1]
        actent <- nl$upperl[acsel + 1]
        cntr <- cntr + 1
        if(action == "S") {
            if(printing == TRUE)
                cat("There are ", cntr, 
                  " reconstruction steps\n")
            return(invisible(list(indexlist = acsellist, rvector = 
                rvector)))
        }
        else if(action == "L")
            acsel <- 2 * acsel
        else {
            acsel <- 2 * acsel + 1
            rvector <- rvector + power
        }
        power <- power * 2
        if(printing == TRUE) {
            cat("Level : ", i, " Action is ", action)
            cat(" (getpacket Index: ", acsel, ")\n")
        }
        acsellist <- c(acsellist, acsel)
    }
    if(printing == TRUE)
        cat("There are ", cntr, " reconstruction steps\n")
    invisible(list(indexlist = acsellist, rvector = rvector))
}
"print.nvwp"<-
function(x, printing = TRUE, ...)
{
    nlev <- nlevelsWT(x)
    pkt <- NULL
    level <- NULL
    decompose <- x$node.list[[nlev]]$upperctrl
    if(decompose == "B") {
        parent.decompose <- 0
        for(i in nlev:1) {
            child.lev <- i - 1
            child.decompose <- sort(c(2 * parent.decompose, 2 * 
                parent.decompose + 1))
            if(child.lev == 0)
                ctrl <- rep("T", 2^nlev)
            else ctrl <- x$node.list[[child.lev]]$upperctrl
            for(j in 1:length(child.decompose)) {
                if(ctrl[child.decompose[j] + 1] == "T") {
                  level <- c(level, child.lev)
                  pkt <- c(pkt, child.decompose[j])
                  if(printing == TRUE)
                    cat("Level: ", child.lev, " Packet: ", 
                      child.decompose[j], "\n")
                }
            }
            if(child.lev != 0) {
                ctrl <- ctrl[child.decompose + 1]
                sv <- ctrl == "B"
                parent.decompose <- child.decompose[sv]
            }
	if (length(parent.decompose)==0)
		break
        }
    }
    else {
        level <- nlev
        pkt <- 0
        if(printing == TRUE) {
            cat("Original data is best packet!\n")
        }
    }
    invisible(list(level = level, pkt = pkt))
}
"print.w2d"<-
function(x, ...)
{
    cat("w2d class object.\n")
    cat("A composite object containing the components\n")
    cat("\t")
    print(names(x))
    cat("Number of levels: ", nlevelsWT(x), "\n")
    cat("Number of data points: ", nrow(x$m), "\n")
    cat("Number of bases: ", ncol(x$m), "\n")
    cat("Groups vector: ")
    print(x$k)
}
"print.wd"<-
function(x, ...)
{
    if(IsEarly(x)) {
        ConvertMessage()
        stop()
    }
    cat("Class 'wd' : Discrete Wavelet Transform Object:\n")
    cat("       ~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    if(x$bc == "interval")
        cat("$transformed.vector is a LONG coefficient vector!\n")
    else cat("$C and $D are LONG coefficient vectors\n")
    cat("\nCreated on :", x$date, "\n")
    cat("Type of decomposition: ", x$type, "\n")
    cat("\nsummary(.):\n----------\n")
    summary.wd(x)
}
"print.wd3D"<-
function(x, ...)
{
    if(IsEarly(x)) {
        ConvertMessage()
        stop()
    }
    cat("Class 'wd3d' : 3D DWT Object:\n")
    cat("       ~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("$ a is the wavelet coefficient array\n")
    cat("Dimension of a is ")
    print(dim(x$a))
    cat("\nCreated on :", x$date, "\n")
    cat("\nsummary(.):\n----------\n")
    summary.wd3D(x)
}
"print.wp"<-
function(x, ...)
{
    if(IsEarly(x)) {
        ConvertMessage()
        stop()
    }
    cat("Class 'wp' : Wavelet Packet Object:\n")
    cat("       ~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("$wp is the wavelet packet matrix\n")
    cat("\nCreated on :", x$date, "\n")
    cat("\nsummary(.):\n----------\n")
    summary.wp(x)
}
"print.wpst"<-
function(x, ...)
{
    if(IsEarly(x)) {
        ConvertMessage()
        stop()
    }
    cat("Class 'wpst' : Nondecimated Wavelet Packet Transform Object:\n")
    cat("       ~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("$wpst is a coefficient vector\n")
    cat("\nCreated on :", x$date[1], "\n")
    cat("\nsummary(.):\n----------\n")
    summary.wpst(x)
}
"print.wpstCL"<-
function(x, ...)
{
    cat("wpstCL class object\n")
    cat("Results of applying discriminator to time series\n")
    cat("Components: ", names(x), "\n")
}
"print.wpstDO"<-
function(x, ...)
{
    cat("Nondecimated wavelet packet discrimination object\n")
    cat("Composite object containing components:")
    print(names(x))
    cat("Fisher's discrimination: done\n")
    cat("BP component has the following information\n")
    print(x$BP)
}
"print.wst"<-
function(x, ...)
{
    if(IsEarly(x)) {
        ConvertMessage()
        stop()
    }
    cat("Class 'wst' : Packet-ordered Nondecimated Wavelet Transform Object:\n")
    cat("       ~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("$wp and $Carray are the coefficient matrices\n")
    cat("\nCreated on :", x$date[1], "\n")
    cat("\nsummary(.):\n----------\n")
    summary.wst(x)
}
"print.wst2D"<-
function(x, ...)
{
    cat("Class 'wst2D' : 2D Packet-ordered Nondecimated Wavelet Transform Object:\n")
    cat("       ~~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("$wst2D is the coefficient array\n")
    cat("\nCreated on :", x$date[1], "\n")
    cat("\nsummary(.):\n----------\n")
    summary.wst2D(x)
}
"putC"<-
function(...)
UseMethod("putC")
"putC.mwd"<-
function(mwd, level, M, boundary = FALSE, index = FALSE, ...)
{
#
#putC.mwd,  changes the C coefficients at the given level.
#Tim Downie
#last update May 1996
#
    if(is.null(class(mwd))) stop("mwd is not class mwd object")
    if(class(mwd) != "mwd")
        stop("mwd is not class mwd object")
    if(level < 0)
        stop("level too small")
    else if(level > nlevelsWT(mwd))
        stop("level too big")
    flc <- mwd$fl.dbase$first.last.c[level + 1,  ]
    if(boundary == FALSE) {
        if(mwd$type == "wavelet")
            n <- 2^level
        else n <- 2^nlevelsWT(mwd)
        i1 <- flc[3] + 1 - flc[1]
        i2 <- flc[3] + n - flc[1]
    }
    else {
        n <- flc[2] - flc[1] + 1
        i1 <- flc[3] + 1
        i2 <- flc[3] + n
    }
    if(index == FALSE) {
        if(length(M) != mwd$filter$npsi * n)
            stop("The length of M is wrong")
        mwd$C[, i1:i2] <- M
        return(mwd)
    }
    else return(list(ix1 = i1, ix2 = i2))
}
"putC.wd"<-
function(wd, level, v, boundary = FALSE, index = FALSE, ...)
{
    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
    if(is.null(class(wd)))
        stop("wd is not class wd object")
    if(class(wd) != "wd")
        stop("wd is not class wd object")
    if(level < 0)
        stop("level should be zero or larger")
    else if(level > nlevelsWT(wd))
        stop(paste("Level should be less than or equal to ", nlevelsWT(wd
            )))
    if(wd$bc == "interval") {
        if(level != wd$current.scale)
            stop(paste(
                "Requested wd object was decomposed to level ", 
                wd$current.scale, 
                " and so for \"wavelets on the interval\" object\ns I can only alter this level for the scaling function coefficients\n"
                ))
        first.level <- wd$fl.dbase$first.last.c[1]
        last.level <- wd$fl.dbase$first.last.c[2]
        offset.level <- wd$fl.dbase$first.last.c[3]
        n <- last.level - first.level + 1
        if(length(v) != n)
            stop(paste(
                "I think the length of \"v\" is wrong. I think it should be of length ",
                n))
        wd$transformed.vector[(offset.level + 1 - first.level):(
            offset.level + n - first.level)] <- v
        return(wd)
    }
    flc <- wd$fl.dbase$first.last.c[level + 1,  ]
    if(boundary == FALSE) {
        if(wd$type == "wavelet")
            n <- 2^level
        else n <- 2^nlevelsWT(wd)
        i1 <- flc[3] + 1 - flc[1]
        i2 <- flc[3] + n - flc[1]
    }
    else {
        n <- flc[2] - flc[1] + 1
        i1 <- flc[3] + 1
        i2 <- flc[3] + n
    }
    if(length(v) != n)
        stop(paste("I think the length of \"v\" is wrong. I think it should be of length ",
            n))
    wd$C[i1:i2] <- v
    if(index == FALSE)
        return(wd)
    else return(list(ix1 = i1, ix2 = i2))
}
"putC.wp"<-
function(wp, ...)
{
    stop("A wavelet packet object does not have ``levels'' of father wavelet coefficients. Use putD to obtain levels of father and mother coefficients"
        )
}
"putC.wst"<-
function(wst, level, value, ...)
{
#
#
# Get all coefficients at a particular level
# First work out how many packets there are at this level
#
    nlevels <- nlevelsWT(wst)
    if(2^nlevels != length(value))
        stop("Input data value of wrong length")
    wst$Carray[level + 1,  ] <- value
    wst
}
"putD"<-
function(...)
UseMethod("putD")
"putD.mwd"<-
function(mwd, level, M, boundary = FALSE, index = FALSE, ...)
{
#
#putD.mwd
#replaces D coefficients at given level with M
#Tim Downie
#last update May 1996
#
#
    if(is.null(class(mwd))) stop("mwd is not class mwd object")
    if(class(mwd) != "mwd")
        stop("mwd is not class mwd object")
    if(level < 0)
        stop("level too small")
    else if(level >= nlevelsWT(mwd))
        stop("level too big")
    fld <- mwd$fl.dbase$first.last.d[level + 1,  ]
    if(boundary == FALSE) {
        if(mwd$type == "wavelet")
            n <- 2^level
        else n <- 2^nlevelsWT(mwd)
        i1 <- fld[3] + 1 - fld[1]
        i2 <- fld[3] + n - fld[1]
    }
    else {
        n <- fld[2] - fld[1] + 1
        i1 <- fld[3] + 1
        i2 <- fld[3] + n
    }
    if(index == FALSE) {
        if(length(M) != mwd$filter$npsi * n)
            stop("The length of M is wrong")
        mwd$D[, i1:i2] <- M
        return(mwd)
    }
    else return(list(ix1 = i1, ix2 = i2))
}
"putD.wd"<-
function(wd, level, v, boundary = FALSE, index = FALSE, ...)
{
    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
    if(is.null(class(wd)))
        stop("wd is not class wd object")
    if(class(wd) != "wd")
        stop("wd is not class wd object")
    if(level < 0)
        stop("level too small")
    else if(level > nlevelsWT(wd)- 1)
        stop(paste("Level too big. Maximum level is ", nlevelsWT(wd)- 1))
    if(wd$bc == "interval") {
        level <- level - wd$current.scale
        objname <- deparse(substitute(wd))
        if(level < 0)
            stop(paste("The wd object: ", objname, 
                " was only decomposed down to level: ", wd$
                current.scale, " Try a larger level"))
        if(boundary == TRUE)
            stop("There are no boundary elements in a wavelets on th\ne interval transform!"
                )
    }
    fld <- wd$fl.dbase$first.last.d[level + 1,  ]
    if(boundary == FALSE) {
        if(wd$type == "wavelet")
            n <- 2^level
        else n <- 2^nlevelsWT(wd)
        if(wd$bc == "interval")
            n <- fld[2] - fld[1] + 1
        i1 <- fld[3] + 1 - fld[1]
        i2 <- fld[3] + n - fld[1]
    }
    else {
        n <- fld[2] - fld[1] + 1
        i1 <- fld[3] + 1
        i2 <- fld[3] + n
    }
    if(length(v) != n)
        stop("I think that the length of v is wrong")
    if(wd$bc == "interval")
        wd$transformed.vector[i1:i2] <- v
    else wd$D[i1:i2] <- v
    if(index == FALSE)
        return(wd)
    else return(list(ix1 = i1, ix2 = i2))
}
"putD.wd3D"<-
function(x, v, ...)
{
    truesize <- dim(x$a)[1]
    nlx <- nlevelsWT(x)
    vlev <- v$lev
    va <- v$a
    putDwd3Dcheck(lti = vlev, dima = dim(va), block = v$block, nlx = nlx)
    Iarrayix <- switch(v$block,
        HHH = 0,
        GHH = 1,
        HGH = 2,
        GGH = 3,
        HHG = 4,
        GHG = 5,
        HGG = 6,
        GGG = 7)
    if(Iarrayix == 0 && vlev != 0)
        stop("Can only insert HHH into level 0")
    if(is.null(Iarrayix))
        stop(paste("Unknown block to insert: ", v$block))
    tmp <- .C("putarr",
        Carray = as.double(x$a),
        truesize = as.integer(truesize),
        level = as.integer(vlev),
        Iarrayix = as.integer(Iarrayix),
        Iarray = as.double(va), PACKAGE = "wavethresh")
    x$a <- array(tmp$Carray, dim = dim(x$a))
    x
}
"putD.wp"<-
function(wp, level, value, ...)
{
#
# Insert coefficients "value" into "wp" at resolution "level".
# First work out how many packets there are at this level
#
    nlev <- nlevelsWT(wp)
    if(2^nlev != length(value))
        stop("Input data value of wrong length")
    wp$wp[level + 1,  ] <- value
    wp
}
"putD.wst"<-
function(wst, level, value, ...)
{
#
#
# Get all coefficients at a particular level
# First work out how many packets there are at this level
#
    nlevels <- nlevelsWT(wst)
    if(2^nlevels != length(value))
        stop("Input data value of wrong length")
    wst$wp[level + 1,  ] <- value
    wst
}
"putDwd3Dcheck"<-
function(lti, dima, block, nlx)
{
    if(lti < 0)
        stop(paste("Level cannot be negative for block:", block))
    else if(lti > nlx - 1)
        stop(paste("Maximum level for block: ", block, " is ", nlx - 1)
            )
    if(length(dima) != 3)
        stop(paste(block, "array is not three-dimensional"))
    if(any(dima != dima[1]))
        stop(paste(block, " dimensions are not all the same"))
    arrdimlev <- IsPowerOfTwo(dima[1])
    if(is.na(arrdimlev))
        stop(paste(block, " dimensions are not power of two"))
    if(arrdimlev != lti)
        stop(paste(block, 
            "dimensions will not fit into cube at that level"))
}
"putpacket"<-
function(...)
UseMethod("putpacket")
"putpacket.wp"<-
function(wp, level, index, packet, ...)
{
#   cat("PUTPACKET: Level:", level, " Index:", index, " Pack Length ", 
#       length(packet), "\n")
    if(class(wp) != "wp") stop("wp object is not of class wp")
    if(level > nlevelsWT(wp))
        stop("Not that many levels in wp object")
    unit <- 2^level
    LocalIndex <- unit * index + 1
    if(index > 2^(nlevelsWT(wp)- level) - 1) {
        cat("Index was too high, maximum for this level is ", 2^(wp$
            nlevels - level) - 1, "\n")
        stop("Error occured")
    }
    if(LocalIndex < 0)
        stop("Index must be  non-negative")
    if(length(packet) != unit)
        stop("Packet is not of correct length\n")
    wp$wp[level + 1, (LocalIndex:(LocalIndex + unit - 1))] <- packet
    wp
}
"putpacket.wst"<-
function(wst, level, index, packet, ...)
{
    class(wst) <- "wp"
    l <- putpacket.wp(wst, level = level, index = index, packet = packet)
    class(l) <- "wst"
    l
}
"putpacket.wst2D"<-
function(wst2D, level, index, type = "S", packet, Ccode = TRUE, ...)
{
    cellength <- 2^level
    nlev <- nlevelsWT(wst2D)
    if(!is.matrix(packet))
        stop("packet should be a matrix")
    nr <- nrow(packet)
    nc <- ncol(packet)
    if(nr != nc)
        stop("packet should be a square matrix")
    else if(nr != cellength)
        stop(paste("packet matrix should be square of dimension ", 
            cellength, " if you're inserting at level ", level, 
            " not ", nr))
    if(level > nlev - 1)
        stop(paste("Maximum level is ", nlev - 1, " you supplied ", 
            level))
    else if(level < 0)
        stop(paste("Minimum level is 0 you supplied ", level))
    if(type != "S" && type != "H" && type != "V" && type != "D")
        stop("Type must be one of S, H, V or D")
    if(nchar(index) != nlev - level)
        stop(paste("Index must be ", nlev - level, 
            " characters long for level ", level))
    for(i in 1:nchar(index)) {
        s1 <- substring(index, i, i)
        if(s1 != "0" && s1 != "1" && s1 != "2" && s1 != "3")
            stop(paste("Character ", i, 
                " in index is not a 0, 1, 2 or 3. It is ", s1))
    }
    if(Ccode == TRUE) {
        ntype <- switch(type,
            S = 0,
            H = 1,
            V = 2,
            D = 3)
        amdim <- dim(wst2D$wst2D)
        ans <- .C("putpacketwst2D",
            am = as.double(wst2D$wst2D),
            d1 = as.integer(amdim[1]),
            d12 = as.integer(amdim[1] * amdim[2]),
            maxlevel = as.integer(nlev - 1),
            level = as.integer(level),
            index = as.integer(index),
            ntype = as.integer(ntype),
            packet = as.double(packet),
            sl = as.integer(nr), PACKAGE = "wavethresh")
        wst2D$wst2D <- array(ans$am, dim = amdim)
    }
    else {
        x <- y <- 0
        ans <- .C("ixtoco",
            level = as.integer(level),
            maxlevel = as.integer(nlev - 1),
            index = as.integer(index),
            x = as.integer(x),
            y = as.integer(y), PACKAGE = "wavethresh")
        tmpx <- switch(type,
            S = 0,
            H = 0,
            V = cellength,
            D = cellength)
        tmpy <- switch(type,
            S = 0,
            H = cellength,
            V = 0,
            D = cellength)
        x <- ans$x + tmpx + 1
        y <- ans$y + tmpy + 1
        cat("x ", x, "y: ", y, "x+cellength-1 ", x + cellength - 1, 
            "y+cellength-1", y + cellength - 1, "\n")
        wst2D$wst2D[level + 1, x:(x + cellength - 1), y:(y + cellength - 
            1)] <- packet
    }
    wst2D
}
"rcov"<-
function(x)
{
#
#rcov
#
#computes a robust correlation matrix of x
# x must be a matrix with the columns as observations
#which is the opposite to the S function var (don't get confused!)
#Method comes from Huber's "Robust Statistics"
#
    if(!is.matrix(x)) stop("x must be a matrix")
    m <- dim(x)[1]
    n <- dim(x)[2]
    b1 <- b2 <- b3 <- 0
    a <- rep(0, m)
    sigma <- matrix(rep(0, m^2), nrow = m)
    for(i in 1:m) {
        a[i] <- 1/mad(x[i,  ])
        sigma[i, i] <- 1/a[i]^2
    }
    if(m > 1) {
        for(i in 2:m)
            for(j in 1:(i - 1)) {
                b1 <- mad(a[i] * x[i,  ] + a[j] * x[j,  ])^2
                b2 <- mad(a[i] * x[i,  ] - a[j] * x[j,  ])^2
                b3 <- mad(a[j] * x[j,  ] - a[i] * x[i,  ])^2
                sigma[i, j] <- (b1 - b2)/((b1 + b2) * a[i] * a[
                  j])
                sigma[j, i] <- (b1 - b3)/((b1 + b3) * a[i] * a[
                  j])
            }
    }
    return(sigma)
}
"rfft"<-
function(x)
{
# given a vector x computes the real continuous fourier transform of
#  x;  ie regards x as points on a periodic function on [0,1] starting at
#  0  and finding the coefficients of the functions 1, sqrt(2)cos 2 pi t, 
#  sqrt(2) sin 2 pi t, etc that give an expansion of the interpolant of
# x    The number of terms in the expansion is the length of x.
# If x is of even length, the last 
#  coefficient will be that of a cosine term with no matching sine.
#
    nx <- length(x)
    z <- fft(x)
    z1 <- sqrt(2) * z[2:(1 + floor(nx/2))]
    rz <- c(Re(z)[1], as.vector(rbind(Re(z1),  - Im(z1))))/nx
    return(rz[1:nx])
}
"rfftinv"<-
function(rz, n = length(rz))
{
#  Inverts the following transform----
# given a vector rz computes the inverse real continuous fourier transform of
#  rz;  ie regards rz as the coefficients of the expansion of a 
#  periodic function f in terms of the functions 
#   1, sqrt(2)cos 2 pi t,   sqrt(2) sin 2 pi t, etc .  
#   The output of the function is f evaluated
# at a regular grid of n points, starting at 0. 
#   If n is not specified it is taken to be the length of rz;
#   the results are unpredictable if n < length(rz).
#
    nz <- length(rz)
    z <- complex(n)
    nz1 <- floor(nz/2)
    nz2 <- ceiling(nz/2) - 1
    z[1] <- rz[1] + (0i)
    z[2:(nz1 + 1)] <- (1/sqrt(2)) * rz[seq(from = 2, by = 2, length = nz1)]
    z[2:(nz2 + 1)] <- z[2:(nz2 + 1)] - (1i) * (1/sqrt(2)) * rz[seq(from = 3,
        by = 2, length = nz2)]
    z[n:(n + 1 - nz1)] <- Conj(z[2:(nz1 + 1)])
    x <- Re(fft(z, inverse = TRUE))
    return(x)
}
"rfftwt"<-
function(xrfft, wt)
{
#    weight the real fourier series xrfft of even length
#     by a weight sequence wt
#    The first term of xrfft is left alone, and the weights are
#    then applied to pairs of terms in xrfft.
#       wt is of length half n .
    xsrfft <- xrfft * c(1, rep(wt, c(rep(2, length(wt) - 1), 1)))
    return(xsrfft)
}
"rm.det"<-
function(wd.int.obj)
{
    len <- length(wd.int.obj$transformed.vector)
    n <- len
    maxscale <- log(len, 2)
    minscale <- wd.int.obj$current.scale
    for(i in c(maxscale:(minscale + 1)))
        n <- n/2
    for(i in c((n + 1):len))
        wd.int.obj$transformed.vector[i] <- 0
    return(wd.int.obj)
}
"rmget"<-
function(requestJ, filter.number, family)
{
    ps <- paste("rm.*.", filter.number, ".", family, sep = "")
    cand <- objects(envir = WTEnv, pattern = ps)
    if(length(cand) == 0)
        return(NULL)
    cand <- substring(cand, first = 4)
    candfd <- firstdot(cand)
    cand <- as.numeric(substring(cand, first = 1, last = candfd - 1))
    cand <- cand[cand >= requestJ]
    if(length(cand) == 0)
        return(NULL)
    else return(min(cand))
}
"rmname"<-
function(J, filter.number, family)
{
    if(J >= 0)
        stop("J must be a negative integer")
    return(paste("rm.",  - J, ".", filter.number, ".", family, sep = ""))
}
"rotateback"<-
function(v)
{
    lv <- length(v)
    v[c(lv, 1:(lv - 1))]
}
"rsswav"<-
function(noisy, value = 1, filter.number = 10, family = "DaubLeAsymm", 
    thresh.type = "hard", ll = 3)
{
    lo <- length(noisy)
    oodd <- noisy[seq(from = 1, by = 2, length = lo/2)]
    oeven <- noisy[seq(from = 2, by = 2, length = lo/2)]    #
#
#   Do decomposition of odd
#
    oddwd <- wd(oodd, filter.number = filter.number, family = family)
    oddwdt <- threshold(oddwd, policy = "manual", value = value, type = 
        thresh.type, lev = ll:(nlevelsWT(oddwd)- 1))
    oddwr <- wr(oddwdt) #
# Interpolate evens
#
    eint <- (c(oeven[1], oeven) + c(oeven, oeven[length(oeven)]))/2
    eint <- eint[1:(length(eint) - 1)]
    ssq1 <- ssq(eint, oddwr)    #
#   ts.plot(oddwr, main = paste("Odd plot, ssq=", ssq1)) #
#   Now do decomposition of even
#
    evenwd <- wd(oeven, filter.number = filter.number, family = family)
    evenwdt <- threshold(evenwd, policy = "manual", value = value, type = 
        thresh.type, lev = ll:(nlevelsWT(evenwd)- 1))
    evenwr <- wr(evenwdt)   #
#
#   Inerpolate odds
#
    oint <- (c(oodd[1], oodd) + c(oodd, oodd[length(oodd)]))/2
    oint <- oint[1:(length(oint) - 1)]
    ssq2 <- ssq(oint, evenwr)   
    #   ts.plot(evenwr, main = paste("Even plot, ssq=", ssq2))
    answd <- wd(noisy, filter.number = filter.number, family = family)
    ll <- list(ssq = (ssq1 + ssq2)/2, df = dof(threshold(answd, policy = 
        "manual", value = value, type = thresh.type, lev = ll:(answd$
        nlevels - 1))))
    return(ll)
}
"simchirp"<-
function(n = 1024)
{
    x <- 1.0000000000000001e-05 + seq(from = -1, to = 1, length = n + 1)[1:
        n]
    y <- sin(pi/x)
    list(x = x, y = y)
}
"ssq"<-
function(u, v)
{
    sum((u - v)^2)
}
"summary.imwd"<-
function(object, ...)

{
#
#
#       Check class of imwd
#
    ctmp <- class(object)
    if(is.null(ctmp))
        stop("imwd has no class")
    else if(ctmp != "imwd")
        stop("imwd is not of class imwd")
    first.last.c <- object$fl.dbase$first.last.c
    pix <- first.last.c[nlevelsWT(object)+ 1, 2] - first.last.c[nlevelsWT(object)+ 
        1, 1] + 1
    cat("UNcompressed image wavelet decomposition structure\n")
    cat("Levels: ", nlevelsWT(object), "\n")
    cat("Original image was", pix, "x", pix, " pixels.\n")
    cat("Filter was: ", object$filter$name, "\n")
    cat("Boundary handling: ", object$bc, "\n")
}
"summary.imwdc"<-
function(object, ...)
{
#
#
#       Check class of imwdc
#
    ctmp <- class(object)
    if(is.null(ctmp))
        stop("imwdc has no class")
    else if(ctmp != "imwdc")
        stop("imwdc is not of class imwdc")
    first.last.c <- object$fl.dbase$first.last.c
    pix <- first.last.c[nlevelsWT(object)+ 1, 2] - first.last.c[nlevelsWT(object)+ 
        1, 1] + 1
    cat("Compressed image wavelet decomposition structure\n")
    cat("Levels: ", nlevelsWT(object), "\n")
    cat("Original image was", pix, "x", pix, " pixels.\n")
    cat("Filter was: ", object$filter$name, "\n")
    cat("Boundary handling: ", object$bc, "\n")
}
"summary.mwd"<-
function(object, ...)
{
    ctmp <- class(object, ...)
    if(is.null(ctmp))
        stop("Input must have class mwd")
    else if(ctmp != "mwd")
        stop("Input must have class mwd")
    cat("Length of original: ", object$ndata, "\n")
    cat("Levels: ", nlevelsWT(object), "\n")
    cat("Filter was: ", object$filter$name, "\n")
    cat("Scaling fns: ", object$filter$nphi, "\n")
    cat("Wavelet fns: ", object$filter$npsi, "\n")
    cat("Prefilter: ", object$prefilter, "\n")
    cat("Scaling factor: ", object$filter$ndecim, "\n")
    cat("Boundary handling: ", object$bc, "\n")
    cat("Transform type: ", object$type, "\n")
    cat("Date: ", object$date, "\n")
}
"summary.wd"<-
function(object, ...)
{
    if(IsEarly(object)) {
        ConvertMessage()
        stop()
    }
    if(object$bc != "interval")
        pix <- length(accessC(object))
    else pix <- 2^nlevelsWT(object)
    cat("Levels: ", nlevelsWT(object), "\n")
    cat("Length of original: ", pix, "\n")
    cat("Filter was: ", object$filter$name, "\n")
    cat("Boundary handling: ", object$bc, "\n")
    if(object$bc == "interval")
        if(object$preconditioned == TRUE)
            cat("Preconditioning is ON\n")
        else cat("Preconditioning is OFF\n")
    cat("Transform type: ", object$type, "\n")
    cat("Date: ", object$date, "\n")
}
"summary.wd3D"<-
function(object, ...)
{
    if(IsEarly(object)) {
        ConvertMessage()
        stop()
    }
    cat("Levels: ", nlevelsWT(object), "\n")
    cat("Filter number was: ", object$filter.number, "\n")
    cat("Filter family was: ", object$family, "\n")
    cat("Date: ", object$date, "\n")
}
"summary.wp"<-
function(object, ...)
{
    if(IsEarly(object)) {
        ConvertMessage()
        stop()
    }
    wpdim <- dim(object$wp)
    cat("Levels: ", nlevelsWT(object), "\n")
    cat("Length of original: ", wpdim[2], "\n")
    cat("Filter was: ", object$filter$name, "\n")
}
"summary.wpst"<-
function(object, ...)
{
    if(IsEarly(object)) {
        ConvertMessage()
        stop()
    }
    pix <- 2^nlevelsWT(object)
    cat("Levels: ", nlevelsWT(object), "\n")
    cat("Length of original: ", pix, "\n")
    cat("Filter was: ", object$filter$name, "\n")
    cat("Date: ", object$date[1], "\n")
    if(length(object$date) != 1)
        cat("This object has been modified. Use \"Whistory\" to find out what's happened\n"
            )
}
"summary.wst"<-
function(object, ...)
{
    if(IsEarly(object)) {
        ConvertMessage()
        stop()
    }
    pix <- 2^nlevelsWT(object)
    cat("Levels: ", nlevelsWT(object), "\n")
    cat("Length of original: ", pix, "\n")
    cat("Filter was: ", object$filter$name, "\n")
    cat("Date: ", object$date[1], "\n")
    if(length(object$date) != 1)
        cat("This object has been modified. Use \"Whistory\" to find out what's happened\n"
            )
}
"summary.wst2D"<-
function(object, ...)
{
    nlev <- nlevelsWT(object)
    cat("Levels: ", nlev, "\n")
    cat("Length of original: ", 2^nlev, "x", 2^nlev, "\n")
    cat("Filter was: ", object$filter$name, "\n")
    cat("Date: ", object$date[1], "\n")
    if(length(object$date) != 1)
        cat("This object has been modified. Use \"Whistory\" to find out what's happened\n"
            )
}
"support"<-
function(filter.number = 10, family = "DaubLeAsymm", m = 0, n = 0)
{
    m <- m + 1
    if(family == "DaubExPhase") {
        a <-  - (filter.number - 1)
        b <- filter.number
        lh <- 2^( + m) * (a + n)
        rh <- 2^( + m) * (b + n)
        return(list(lh = lh, rh = rh, psi.lh =  - (filter.number - 1), 
            psi.rh = filter.number, phi.lh = 0, phi.rh = 2 * 
            filter.number - 1))
    }
    else if(family == "DaubLeAsymm") {
        a <-  - (filter.number - 1)
        b <- filter.number
        lh <- 2^( + m) * (a + n)
        rh <- 2^( + m) * (b + n)
        return(list(lh = lh, rh = rh, psi.lh =  - (filter.number - 1), 
            psi.rh = filter.number, phi.lh = 0, phi.rh = 2 * 
            filter.number - 1))
    }
    else {
        stop(paste("Family: ", family, " not supported for support!\n")
            )
    }
}
"sure"<-
function(x)
{
#
# The SURE function of Donoho and Johnstone
# Finds the minimum
#
    x <- abs(x)
    d <- length(x)
    y <- sort(x)    #
#
#       Form cumulative sum
#
    cy <- cumsum(y^2)
    cy <- c(0, cy[1:(length(cy) - 1)])  #
#
#       Now the answer
#
    ans <- d - 2 * 1:d + cy + d:1 * y^2 #       cat("ans is\n")
#       print(ans)
    m <- min(ans)
    index <- (1:length(ans))[m == ans]
    return(y[index])
}
"threshold"<-
function(...)
UseMethod("threshold")
"threshold.imwd"<-
function(imwd, levels = 3:(nlevelsWT(imwd)- 1), type = "hard", policy = 
    "universal", by.level = FALSE, value = 0, dev = var, verbose = FALSE, 
    return.threshold = FALSE, compression = TRUE, Q = 0.050000000000000003, ...)
{
#
#
#   Check class of imwd
#
    if(verbose == TRUE) cat("Argument checking\n")
    ctmp <- class(imwd)
    if(is.null(ctmp))
        stop("imwd has no class")
    else if(ctmp != "imwd")
        stop("imwd is not of class imwd")
    if(policy != "universal" && policy != "manual" && policy != 
        "probability" && policy != "fdr")
        stop("Only policys are universal, manual, fdr and probability at present"
            )
    if(type != "hard" && type != "soft")
        stop("Only hard or soft thresholding at  present")
    r <- range(levels)
    if(r[1] < 0)
        stop("levels out of range, level too small")
    if(r[2] > nlevelsWT(imwd)- 1)
        stop("levels out of range, level too big")
    if(r[1] > nlevelsWT(imwd)- 1) {
        warning("no thresholding done")
        return(imwd)
    }
    if(r[2] < 0) {
        warning("no thresholding done")
        return(imwd)
    }
    nthresh <- length(levels)
    d <- NULL
    n <- 2^(2 * nlevelsWT(imwd))  #
#       Decide which policy to adopt
#               The next if-else construction should define a vector called
#               "thresh" that contains the threshold value for each level
#               in "levels". This may be the same threshold value
#               a global threshold.
#
    if(policy == "universal") {
        if(verbose == TRUE)
            cat("Universal policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh) {
                d <- c(d, imwd[[lt.to.name(levels[i], "CD")]], 
                  imwd[[lt.to.name(levels[i], "DC")]], imwd[[
                  lt.to.name(levels[i], "DD")]])
            }
            noise.level <- sqrt(dev(d))
            thresh <- sqrt(2 * log(n)) * noise.level
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh, "\n")
            thresh <- rep(thresh, length = nthresh)
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for(i in 1:nthresh) {
                d <- c(imwd[[lt.to.name(levels[i], "CD")]], 
                  imwd[[lt.to.name(levels[i], "DC")]], imwd[[
                  lt.to.name(levels[i], "DD")]])
                noise.level <- sqrt(dev(d))
                thresh[i] <- sqrt(2 * log(n)) * noise.level
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
    }
    else if(policy == "manual") {
        if(verbose == TRUE)
            cat("Manual policy...\n")
        thresh <- rep(value, length = nthresh)
        if(length(value) != 1 && length(value) != nthresh)
            warning("your threshold is not the same length as number of levels"
                )
    }
    else if(policy == "fdr") {
#
#
#               Threshold chosen by FDR-procedure
#
        if(verbose == TRUE) cat("FDR policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh) {
                d <- c(d, imwd[[lt.to.name(levels[i], "CD")]], 
                  imwd[[lt.to.name(levels[i], "DC")]], imwd[[
                  lt.to.name(levels[i], "DD")]])
            }
            if(length(value) != 1)
                stop("Length of value should be 1")
            noise.level <- sqrt(dev(c(imwd[[lt.to.name(levels[
                nthresh], "CD")]], imwd[[lt.to.name(levels[
                nthresh], "DC")]], imwd[[lt.to.name(levels[
                nthresh], "DD")]])))
            minit <- n
            dinit <- d
            thinit <- qnorm(1 - Q/2) * noise.level
            if(log(n, 2) > 15)
                ninit <- 4
            else {
                if(log(n, 2) > 12)
                  ninit <- 3
                else {
                  if(log(n, 2) > 10)
                    ninit <- 2
                  else ninit <- 1
                }
            }
            for(k in seq(1, ninit)) {
                dinit1 <- dinit[abs(dinit) >= thinit]
                minit <- length(dinit1)
                if(minit == 0)
                  thresh <- max(abs(d)) * 1.0001
                else {
                  thinit <- qnorm(1 - (Q * minit)/(2 * n)) * 
                    noise.level
                  minit1 <- length(dinit1[abs(dinit1) >= thinit
                    ])
                  if(minit1 == minit || minit1 == 0)
                    break
                  dinit <- dinit1
                }
            }
            if(noise.level > 0) {
                m <- length(d)
                minit <- length(dinit)
                p <- (2 - 2 * pnorm(abs(dinit)/noise.level))
                index <- order(p)
                j <- seq(1, minit)
                m0 <- max(j[p[index] <= (Q * j)/m])
                if(m0 != "NA" && m0 < minit)
                  thresh <- abs(dinit[index[m0]])
                else {
                  if(m0 == "NA")
                    thresh <- max(abs(dinit)) * 1.0001
                  else thresh <- 0
                }
            }
            else thresh <- 0
            thresh <- rep(thresh, length = nthresh)
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh[1], "\n", 
                  "sigma is: ", noise.level, "\n")
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for(i in 1:nthresh) {
                d <- c(imwd[[lt.to.name(levels[i], "CD")]], 
                  imwd[[lt.to.name(levels[i], "DC")]], imwd[[
                  lt.to.name(levels[i], "DD")]])
                m <- length(d)
                noise.level <- sqrt(dev(d))
                thinit <- qnorm(1 - Q/2) * noise.level
                dinit <- d[abs(d) >= thinit]
                minit <- length(dinit)
                if(minit == 0)
                  thresh[i] <- max(abs(d)) * 1.0001
                else {
                  if(noise.level > 0) {
                    p <- (2 - 2 * pnorm(abs(dinit)/noise.level)
                      )
                    index <- order(p)
                    j <- seq(1, minit)
                    m0 <- max(j[p[index] <= (Q * j)/m])
                    if(m0 != "NA" && m0 < minit)
                      thresh[i] <- abs(dinit[index[m0]])
                    else {
                      if(m0 == "NA")
                        thresh[i] <- max(abs(dinit)) * 1.0001
                      else thresh[i] <- 0
                    }
                  }
                  else thresh[i] <- 0
                }
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], "is", 
                    thresh[i], "\n")
            }
        }
    }
    else if(policy == "probability") {
        if(verbose == TRUE)
            cat("Probability policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh) {
                d <- c(d, imwd[[lt.to.name(levels[i], "CD")]], 
                  imwd[[lt.to.name(levels[i], "DC")]], imwd[[
                  lt.to.name(levels[i], "DD")]])
            }
            if(length(value) != 1)
                stop("Length of value should be 1")
            thresh <- rep(quantile(abs(d), prob = value), length = 
                nthresh)
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh[1], "\n")
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            if(length(value) == 1)
                value <- rep(value, nthresh)
            if(length(value) != nthresh)
                stop("Wrong number of probability values")
            for(i in 1:nthresh) {
                d <- c(imwd[[lt.to.name(levels[i], "CD")]], 
                  imwd[[lt.to.name(levels[i], "DC")]], imwd[[
                  lt.to.name(levels[i], "DD")]])
                thresh[i] <- quantile(abs(d), prob = value[i])
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
    }
    if(return.threshold == TRUE)
        return(thresh)
    for(i in 1:nthresh) {
        dCD <- imwd[[lt.to.name(levels[i], "CD")]]
        dDC <- imwd[[lt.to.name(levels[i], "DC")]]
        dDD <- imwd[[lt.to.name(levels[i], "DD")]]
        if(type == "hard") {
            dCD[abs(dCD) <= thresh[i]] <- 0
            dDC[abs(dDC) <= thresh[i]] <- 0
            dDD[abs(dDD) <= thresh[i]] <- 0
            if(verbose == TRUE) {
                cat("Level: ", levels[i], " there are ", sum(
                  dCD == 0), ":", sum(dDC == 0), ":", sum(dDD == 
                  0), " zeroes and: ")
                cat(sum(dCD != 0), ":", sum(dDC != 0), ":", sum(
                  dDD != 0), " nonzeroes\n")
            }
        }
        else if(type == "soft") {
            dCD <- sign(dCD) * (abs(dCD) - thresh[i]) * (abs(dCD) > 
                thresh[i])
            dDC <- sign(dDC) * (abs(dDC) - thresh[i]) * (abs(dDC) > 
                thresh[i])
            dDD <- sign(dDD) * (abs(dDD) - thresh[i]) * (abs(dDD) > 
                thresh[i])
            if(verbose == TRUE) {
                cat("Level: ", levels[i], " there are ", sum(
                  dCD == 0), ":", sum(dDC == 0), ":", sum(dDD == 
                  0), " zeroes and: ")
                cat(sum(dCD != 0), ":", sum(dDC != 0), ":", sum(
                  dDD != 0), " nonzeroes\n")
            }
        }
        imwd[[lt.to.name(levels[i], "CD")]] <- dCD
        imwd[[lt.to.name(levels[i], "DC")]] <- dDC
        imwd[[lt.to.name(levels[i], "DD")]] <- dDD
    }
    if(compression == TRUE)
        return(compress(imwd, verbose = verbose))
    else return(imwd)
}
"threshold.imwdc"<-
function(imwdc, verbose = FALSE, ...)
{
    warning("You are probably thresholding an already thresholded object")
    imwd <- uncompress(imwdc, verbose = verbose)
    return(threshold(imwd, verbose = TRUE, ...))
}
"threshold.irregwd"<-
function(irregwd, levels = 3:(nlevelsWT(wd)- 1), type = "hard", policy = 
    "universal", by.level = FALSE, value = 0, dev = var, boundary = FALSE, verbose
     = FALSE, return.threshold = FALSE, force.sure = FALSE, cvtol = 0.01, Q = 
    0.050000000000000003, alpha = 0.050000000000000003, ...)
{
    if(verbose == TRUE)
        cat("threshold.irregwd:\n")
    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
#
#   Check class of wd
#
    if(verbose == TRUE)
        cat("Argument checking\n")
    ctmp <- class(irregwd)
    if(is.null(ctmp))
        stop("irregwd has no class")
    else if(ctmp != "irregwd")
        stop("irregwd is not of class irregwd")
    wd <- irregwd
    class(wd) <- "wd"
    if(policy != "universal" && policy != "manual" && policy != 
        "probability" && policy != "sure" && policy != "mannum" && 
        policy != "cv" && policy != "fdr" && policy != "op1" && policy != 
        "op2" && policy != "LSuniversal")
        stop("Only policys are universal, manual, mannum, sure, LSuniversal, cv, op1, op2 and probability at present"
            )
    if(type != "hard" && type != "soft")
        stop("Only hard or soft thresholding at  present")
    r <- range(levels)
    if(r[1] < 0)
        stop("levels out of range, level too small")
    if(r[2] > nlevelsWT(wd)- 1)
        stop("levels out of range, level too big")
    if(r[1] > nlevelsWT(wd)- 1) {
        warning("no thresholding done")
        return(wd)
    }
    if(r[2] < 0) {
        warning("no thresholding done")
        return(wd)
    }
    n <- 2^nlevelsWT(wd)
    nthresh <- length(levels)   #
# Estimate sigma
    if(by.level == FALSE) {
        d <- NULL
        ccc <- NULL
        for(i in 1:nthresh) {
            d <- c(d, accessD(wd, level = levels[i], boundary = 
                boundary))
            ccc <- c(ccc, accessc(irregwd, level = levels[i], 
                boundary = boundary))
        }
        ind <- (1:length(d))[abs(ccc) > 1.0000000000000001e-05]
        sigma <- sqrt(dev(d[ind]/sqrt(ccc[ind])))
        sigma <- rep(sigma, nthresh)
    }
    else {
        for(i in 1:nthresh) {
            d <- accessD(wd, level = levels[i], boundary = boundary
                )
            ccc <- accessc(irregwd, level = levels[i], boundary = 
                boundary)
            ind <- (1:length(d))[abs(ccc) > 1.0000000000000001e-05]
            sigma[i] <- sqrt(dev(d[ind]/sqrt(ccc[ind])))
        }
    }
    if(verbose == TRUE)
        print(sigma)
    d <- NULL
    ccc <- NULL #
#   Check to see if we're thresholding a complex wavelet transform.
#   We can only do certain things in this case
#
    if(is.complex(wd$D)) {
        stop("Complex transform not implemented")
    }
#
#
#   Decide which policy to adopt
#       The next if-else construction should define a vector called
#       "thresh" that contains the threshold value for each level
#       in "levels". This may be the same threshold value
#       a global threshold.
#
    if(policy == "universal") {
#
#
#       Donoho and Johnstone's universal policy
#
        if(verbose == TRUE) cat("Universal policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh)
                d <- c(d, accessD(wd, level = levels[i], 
                  boundary = boundary))
            nd <- length(d)
            thresh <- sqrt(2 * log(nd))
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh, "\n")
            thresh <- rep(thresh, length = nthresh)
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for(i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = 
                  boundary)
                nd <- length(d)
                thresh[i] <- sqrt(2 * log(nd))
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
        expo <- 1
    }
    else if(policy == "LSuniversal") {
#
#
#       The universal policy modified for local spectral smoothing
#       This should only be used via the LocalSpec function
#
        if(verbose == TRUE) cat("Local spectral universal policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh)
                d <- c(d, accessD(wd, level = levels[i], 
                  boundary = boundary))
            nd <- length(d)
            thresh <- log(nd)
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh, "\n")
            thresh <- rep(thresh, length = nthresh)
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for(i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = 
                  boundary)
                nd <- length(d)
                thresh[i] <- log(nd)
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
        expo <- 1
    }
    else if(policy == "sure") {
        if(type == "hard")
            stop("Can only do soft thresholding with sure policy")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh) {
                d <- c(d, accessD(wd, level = levels[i], 
                  boundary = boundary))
                ccc <- c(ccc, accessc(irregwd, level = levels[i
                  ], boundary = boundary))
            }
            ind <- (1:length(d))[abs(ccc) > 1.0000000000000001e-05]
            nd <- length(ind)
            neta.d <- (log(nd, base = 2)^(3/2))
            sd2 <- (sum((d[ind]/(sigma[1] * ccc)[ind])^2 - 1)/nd)
            if(verbose == TRUE) {
                cat("neta.d is ", neta.d, "\nsd2 is ", sd2, 
                  "\n")
                cat("nd is ", nd, "\n")
                cat("noise.level ", noise.level, "\n")
            }
            if(force.sure == TRUE || sd2 > neta.d/sqrt(nd)) {
                if(verbose == TRUE) {
                  cat("SURE: Using SURE\n")
                }
                thresh <- newsure(sqrt(ccc) * sigma[1], d)
                expo <- 0
            }
            else {
                if(verbose == TRUE)
                  cat("SURE: (sparse) using sqrt 2log n\n")
                thresh <- sqrt(2 * log(nd))
            }
            thresh <- rep(thresh, length = nthresh)
            if(verbose == TRUE)
                cat("Global threshold is ", thresh, "\n")
        }
        else {
#
#
#       By level is true
#
            print("Sure for level- and coefficient-dependenet thresholding is not adapted"
                )
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            collect <- NULL
            for(i in 1:nthresh)
                collect <- c(collect, accessD(wd, level = 
                  levels[i], boundary = boundary))
            noise.level <- sqrt(dev(collect))
            for(i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = 
                  boundary)
                nd <- length(d)
                neta.d <- (log(nd, base = 2)^(3/2))
                sd2 <- (sum((d/noise.level)^2 - 1)/nd)
                if(verbose == TRUE) {
                  cat("neta.d is ", neta.d, "\nsd2 is ", sd2, 
                    "\n")
                  cat("nd is ", nd, "\n")
                  cat("noise.level ", noise.level, "\n")
                }
                if(force.sure == TRUE || sd2 > neta.d/sqrt(nd)) {
                  if(verbose == TRUE) {
                    cat("SURE: Using SURE\n")
                  }
                  thresh[i] <- sure(d/noise.level)
                }
                else {
                  if(verbose == TRUE)
                    cat("SURE: (sparse) using sqrt 2log n\n")
                  thresh[i] <- sqrt(2 * log(nd))
                }
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
    }
    else if(policy == "manual") {
#
#
#       User supplied threshold policy
#
        if(verbose == TRUE) cat("Manual policy\n")
        thresh <- rep(value, length = nthresh)
        expo <- 1
        if(length(value) != 1 && length(value) != nthresh)
            warning("your threshold is not the same length as number of levels"
                )
    }
    else if(policy == "mannum") {
        if(verbose == TRUE) {
            cat("Manual policy using ", value, " of the")
            cat(" largest coefficients\n")
        }
        if(value < 1) {
            stop("Have to select an integer larger than 1 for value"
                )
        }
        else if(value > length(wd$D)) {
            stop(paste("There are only ", length(wd$D), 
                " coefficients, you specified ", value))
        }
        coefs <- wd$D
        scoefs <- sort(abs(coefs))
        scoefs <- min(rev(scoefs)[1:value])
        wd$D[abs(wd$D) < scoefs] <- 0
        return(wd)
    }
    else if(policy == "probability") {
#
#
#       Threshold is quantile based
#
        if(verbose == TRUE) cat("Probability policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh)
                d <- c(d, accessD(wd, level = levels[i], 
                  boundary = boundary))
            if(length(value) != 1)
                stop("Length of value should be 1")
            thresh <- rep(quantile(abs(d), prob = value), length = 
                nthresh)
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh[1], "\n")
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            if(length(value) == 1)
                value <- rep(value, nthresh)
            if(length(value) != nthresh)
                stop("Wrong number of probability values")
            for(i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = 
                  boundary)
                thresh[i] <- quantile(abs(d), prob = value[i])
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
    }
    if(return.threshold == TRUE)
        return(thresh)
    for(i in 1:nthresh) {
        d <- accessD(wd, level = levels[i], boundary = boundary)
        ccc <- accessc(irregwd, level = levels[i], boundary = boundary)
        actthresh <- thresh[i] * (sigma[i] * sqrt(ccc))^expo    
    # is vector
        if(type == "hard") {
            d[abs(d) <= actthresh] <- 0
            if(verbose == TRUE)
                cat("Level: ", levels[i], " there are ", sum(d == 
                  0), " zeroes\n")
        }
        else if(type == "soft") {
            d <- (d * (abs(d) - actthresh) * (abs(d) > actthresh))/
                abs(d)
            d[is.na(d)] <- 0
        }
        wd <- putD(wd, level = levels[i], v = d, boundary = boundary)
    }
    wd
}
"threshold.mwd"<-
function(mwd, levels = 3:(nlevelsWT(mwd)- 1), type = "hard", policy = "universal",
    boundary = FALSE, verbose = FALSE, return.threshold = FALSE, threshold = 0, covtol
     = 1.0000000000000001e-09, robust = TRUE, return.chisq = FALSE, bivariate = TRUE, ...)
{
#threshold.mwd
#thresholds a multiple wavelet object
#Tim Downie
#last updated May 1996
#
#
#   Check arguments
#
    if(verbose == TRUE) cat("threshold.mwd:\n")
    if(verbose == TRUE)
        cat("Argument checking\n")
    ctmp <- class(mwd)
    if(is.null(ctmp))
        stop("mwd has no class")
    else if(ctmp != "mwd")
        stop("mwd is not of class mwd")
    if(policy != "manual" && policy != "universal" && policy != 
        "visushrink")
        stop("Only policies are universal manual and visushrink at present"
            )
    if(type != "hard" && type != "soft")
        stop("Only hard or soft thresholding at present")
    nlevels <- nlevelsWT(mwd)
    npsi <- mwd$filter$npsi
    r <- range(levels)
    if(r[1] < 0)
        stop("levels out of range, level too small")
    if(r[2] > nlevelsWT(mwd)- 1)
        stop("levels out of range, level too big")
    if(r[1] > nlevelsWT(mwd)- 1) {
        warning("no thresholding done, returning input")
        return(mwd)
    }
    if(r[2] < 0) {
        warning("no thresholding done, returning input")
        return(mwd)
    }
    if(policy == "manual" && threshold <= 0) stop(
            "If you want manual thresholding, you must supply\na positive threshold"
            )   #
#
#Apply the a single wavelet policy to multiwavelets   
#so far only universal thresholding 
#visushrink visushrink can be done if using the single policy
#
    if(bivariate == FALSE) {
        if(verbose == TRUE)
            cat("Thresholding multiple wavelets using single wavelet universal thresholding\n"
                )
        noise.level <- rep(0, npsi)
        thresh <- rep(0, npsi)
        ninlev <- rep(0, length(levels))
        if(robust == FALSE)
            dev <- var
        else dev <- mad
        D <- NULL
        for(i in levels) {
            index <- i + 1 - levels[1]
            ninlev[index] <- dim(accessD(mwd, level = i, boundary
                 = boundary))[2]
            D <- matrix(c(D, accessD(mwd, level = i, boundary = 
                boundary)), nrow = npsi)
        }
        nD <- dim(D)[2]
        for(i in 1:npsi) {
            noise.level[i] <- sqrt(dev(D[i,  ]))
            if(policy == "visushrink")
                thresh[i] <- (sqrt(2 * log(nD)) * noise.level[i
                  ])/sqrt(nD)
            else if(policy == "manual")
                thresh[i] <- threshold[i]
            else thresh[i] <- (sqrt(2 * log(nD)) * noise.level[i])
        }
        if(verbose == TRUE) {
            cat("Threshold for each wavelet is: ", thresh, "\n")
            cat("noise levels are : ", noise.level, "\n")
        }
        for(i in 1:npsi) {
            d <- D[i,  ]
            if(type == "hard") {
                d[abs(d) <= thresh[i]] <- 0
            }
            else if(type == "soft") {
                d <- sign(d) * (abs(d) - thresh[i]) * (abs(d) > 
                  thresh[i])
            }
            D[i,  ] <- d
        }
        jump <- 1
        for(i in levels) {
            index <- i + 1 - levels[1]
            mwd <- putD(mwd, level = i, M = D[, jump:(jump + ninlev[
                index] - 1)], boundary = boundary)
            jump <- jump + ninlev[index]
        }
        if(return.threshold == TRUE)
            return(thresh)
        else return(mwd)
    }
#
#
#If we get here then do Multivariate thresholding
# 
    if(policy == "universal" || policy == "manual") {
        n <- 0
        nj <- rep(0, length(levels))
        chisq <- NULL
        chisqkeep <- NULL
        chisqnewkeep <- NULL
        for(i in 1:length(levels)) {
            level <- levels[i]
            d <- accessD(mwd, level = level)
            nj[i] <- dim(d)[2]
            Y <- rep(0, nj[i])  
    # VHAT is the Var/Covar matrix of the data at each level
# estinated using normal estimates or robust estimates
#
            if(robust == FALSE)
                VHAT <- var(t(d))
            if(robust == TRUE) VHAT <- rcov(d)  #
# If the smallest eigen value of VHAT is less than covtol
# we may run into problems when inverting VHAT
# so code chisq as -1 and return the same vector coeff as was input
#
            if(min(abs(eigen(VHAT, only.values = TRUE)$values)) < 
                covtol) {
                warning(paste(
                  "singular variance structure in level ", 
                  level, "this level not thresholded"))
                Y <- rep(-1, nj[i])
            }
            else {
                VINV <- solve(VHAT)
                for(s in 1:npsi)
                  Y <- Y + d[s,  ]^2 * VINV[s, s]
                for(s in 2:npsi)
                  for(t in 1:(s - 1))
                    Y <- Y + 2 * d[s,  ] * d[t,  ] * VINV[s, t]
                n <- n + nj[i]  #
# The above line means that the threshold is caculated using only
# the thresholdable coefficients.
            }
            chisq <- c(chisq, Y)
        }
    }
    if(policy != "manual")
        chithresh <- 2 * log(n)
    else chithresh <- threshold
    if(return.threshold == TRUE) {
        return(chithresh)
    }
    if(return.chisq == TRUE)
        return(chisq)
    lc <- length(chisq)
    dnew <- matrix(rep(0, 2 * lc), nrow = 2)
    d <- NULL
    for(i in 1:length(levels)) {
        d <- matrix(c(d, accessD(mwd, level = levels[i])), nrow = 2)
    }
    if(type == "hard") {
        for(i in 1:lc) {
            keep <- 1 * ((chisq[i] >= chithresh) || (chisq[i] == -1
                ))
            dnew[, i] <- d[, i] * keep
        }
    }
    if(type == "soft") {
        for(i in 1:lc) {
            if(chisq[i] != -1)
                chisqnew <- max(chisq[i] - chithresh, 0)
            if(chisq[i] > 0)
                shrink <- (max(chisq[i] - chithresh, 0))/chisq[
                  i]
            else shrink <- 0
            dnew[, i] <- d[, i] * shrink
        }
    }
    low <- 1
    for(i in 1:length(levels)) {
        mwd <- putD(mwd, level = levels[i], M = dnew[, low:(low - 1 + 
            nj[i])])
        low <- low + nj[i]
    }
    if(verbose == TRUE)
        cat("returning wavelet decomposition\n")
    return(mwd)
}
"threshold.wd"<-
function(wd, levels = 3:(nlevelsWT(wd)- 1), type = "soft", policy = "sure", 
    by.level = FALSE, value = 0, dev = madmad, boundary = FALSE, verbose = FALSE, 
    return.threshold = FALSE, force.sure = FALSE, cvtol = 0.01,
	cvmaxits=500, Q = 
    0.050000000000000003, OP1alpha = 0.050000000000000003, alpha = 0.5, 
    beta = 1, C1 = NA, C2 = NA, C1.start = 100, al.check=TRUE, ...)
{
    if(verbose == TRUE)
        cat("threshold.wd:\n")
    if(IsEarly(wd)) {
        ConvertMessage()
        stop()
    }
#
#   Check class of wd
#
    if(verbose == TRUE)
        cat("Argument checking\n")
    ctmp <- class(wd)
    if(is.null(ctmp))
        stop("wd has no class")
    else if(ctmp != "wd")
        stop("wd is not of class wd")
    if(policy != "universal" && policy != "manual" && policy != 
        "probability" && policy != "sure" && policy != "mannum" && 
        policy != "cv" && policy != "fdr" && policy != "op1" && policy != 
        "op2" && policy != "LSuniversal" && policy != "BayesThresh")
        stop("Only policys are universal, BayesThresh, manual, mannum, sure, LSuniversal, cv, op1, op2 and probability at present"
            )
    if(type != "hard" && type != "soft")
        stop("Only hard or soft thresholding at  present")
    r <- range(levels)
    if(r[1] < 0)
        stop("levels out of range, level too small. Minimum level is 0"
            )
    if(r[2] > nlevelsWT(wd) - 1)
        stop(paste("levels out of range, level too big. Maximum level is",
            nlevelsWT(wd) - 1))
    if(r[1] > nlevelsWT(wd)- 1) {
        warning("no thresholding done")
        return(wd)
    }
    if(r[2] < 0) {
        warning("no thresholding done")
        return(wd)
    }
    if (al.check==TRUE)
	if (all(sort(levels)==levels)==FALSE)
		warning("Entries in levels vector are not ascending. Please check this is what you intend. If so, you can turn this warning off with al.check argument")
    d <- NULL
    n <- 2^nlevelsWT(wd)
    nthresh <- length(levels)   #
#
#   Check to see if we're thresholding a complex wavelet transform.
#   We can only do certain things in this case
#
    if(is.complex(wd$D)) {
		
	stop("Please use cthresh package for complex-valued wavelet shrinkage")
    }
#
#
#   Decide which policy to adopt
#       The next if-else construction should define a vector called
#       "thresh" that contains the threshold value for each level
#       in "levels". This may be the same threshold value
#       a global threshold.
#
    if(policy == "universal") {
#
#
#       Donoho and Johnstone's universal policy
#
        if(verbose == TRUE) cat("Universal policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh)
                d <- c(d, accessD(wd, level = levels[i], 
                  boundary = boundary))
            noise.level <- sqrt(dev(d))
            nd <- length(d)
            thresh <- sqrt(2 * log(nd)) * noise.level
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh, "\n")
            thresh <- rep(thresh, length = nthresh)
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for(i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = 
                  boundary)
                noise.level <- sqrt(dev(d))
                nd <- length(d)
                thresh[i] <- sqrt(2 * log(nd)) * 
                    noise.level
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
    }
    else if(policy == "LSuniversal") {
#
#
#       The universal policy modified for local spectral smoothing
#       This should only be used via the LocalSpec function
#
        if(verbose == TRUE) cat("Local spectral universal policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh)
                d <- c(d, accessD(wd, level = levels[i], 
                  boundary = boundary))
            noise.level <- sqrt(dev(d))
            nd <- length(d)
            thresh <- log(nd) * noise.level
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh, "\n")
            thresh <- rep(thresh, length = nthresh)
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for(i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = 
                  boundary)
                noise.level <- sqrt(dev(d))
                nd <- length(d)
                thresh[i] <- log(nd) * noise.level
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
    }
    else if(policy == "sure") {
        if(type == "hard")
            stop("Can only do soft thresholding with sure policy")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh)
                d <- c(d, accessD(wd, level = levels[i], 
                  boundary = boundary))
            noise.level <- sqrt(dev(d))
            nd <- length(d)
            neta.d <- (log(nd, base = 2)^(3/2))
            sd2 <- (sum((d/noise.level)^2 - 1)/nd)
            if(verbose == TRUE) {
                cat("neta.d is ", neta.d, "\nsd2 is ", sd2, 
                  "\n")
                cat("nd is ", nd, "\n")
                cat("noise.level ", noise.level, "\n")
            }
            if(force.sure == TRUE || sd2 > neta.d/sqrt(nd)) {
                if(verbose == TRUE) {
                  cat("SURE: Using SURE\n")
                }
                thresh <- sure(d/noise.level)
            }
            else {
                if(verbose == TRUE)
                  cat("SURE: (sparse) using sqrt 2log n\n")
                thresh <- sqrt(2 * log(nd))
            }
            thresh <- rep(thresh * noise.level, length = nthresh)
            if(verbose == TRUE)
                cat("Global threshold is ", thresh, "\n")
        }
        else {
#
#
#       By level is true
#
            if(verbose == TRUE) cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            collect <- NULL
            for(i in 1:nthresh)
                collect <- c(collect, accessD(wd, level = 
                  levels[i], boundary = boundary))
            noise.level <- sqrt(dev(collect))
            for(i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = 
                  boundary)
                nd <- length(d)
                neta.d <- (log(nd, base = 2)^(3/2))
                sd2 <- (sum((d/noise.level)^2 - 1)/nd)
                if(verbose == TRUE) {
                  cat("neta.d is ", neta.d, "\nsd2 is ", sd2, 
                    "\n")
                  cat("nd is ", nd, "\n")
                  cat("noise.level ", noise.level, "\n")
                }
                if(force.sure == TRUE || sd2 > neta.d/sqrt(nd)) {
                  if(verbose == TRUE) {
                    cat("SURE: Using SURE\n")
                  }
                  thresh[i] <- sure(d/noise.level)
                }
                else {
                  if(verbose == TRUE)
                    cat("SURE: (sparse) using sqrt 2log n\n")
                  thresh[i] <- sqrt(2 * log(nd))
                }
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
    }
    else if(policy == "BayesThresh") {
#
# Check that all hyperparameters of the prior are non-negative
#
        if(alpha < 0) stop("parameter alpha is negative")
        if(beta < 0)
            stop("parameter beta is negative")
        nthresh <- length(levels)
        nsignal <- rep(0, nthresh)
        noise.level <- sqrt(dev(accessD(wd, level = (nlevelsWT(wd)- 1))))
        v <- 2^( - alpha * levels)
        if(is.na(C1)) {
#
# Estimation of C1 and C2 via universal threshodling
#
            if(C1.start < 0) stop("C1.start is negative")
            universal <- threshold(wd, policy = "universal", type
                 = "hard", dev = dev, by.level = FALSE, levels = 
                levels)
            sum2 <- rep(0, nthresh)
            for(i in 1:nthresh) {
                dun <- accessD(universal, level = levels[i])
                nsignal[i] <- sum(abs(dun) > 10^-10)
                if(nsignal[i] > 0)
                  sum2[i] <- sum(dun[abs(dun) > 0]^2)
            }
            if(sum(nsignal) == 0) {
                wd <- nullevels(wd, levelstonu = levels)
                if(verbose == TRUE)
                  cat(
                    "hyperparameters of the prior are: alpha = ",
                    alpha, "C1 = 0", "beta = ", beta, 
                    "C2 = 0\n")
                return(wd)
            }
            else {
		 fntoopt <- function(C, nsignal, noise.level, wd, sum2, v)				{
			ans<- nsignal * (log(noise.level^2 + C^2 * 
			  v) - 2 * log(pnorm(( - noise.level * sqrt(2 * 
			  log(2^nlevelsWT(wd))))/sqrt(noise.level^2 + C^2 * 
			  v)))) + sum2/(noise.level^2 + C^2 * v)
			sum(ans)
			
			}

		C1 <- optimize(f=fntoopt, interval=c(0, 50*sqrt(C1.start)), 
			nsignal=nsignal, noise.level=noise.level, wd=wd, sum2=sum2, v=v)$minimum^2	
		}
	}
        if(C1 < 0)
            stop("parameter C1 is negative")
        tau2 <- C1 * v
        if(is.na(C2)) {
            p <- 2 * pnorm(( - noise.level * sqrt(2 * log(2^wd$
                nlevels)))/sqrt(noise.level^2 + tau2))
            if(beta == 1)
                C2 <- sum(nsignal/p)/nlevelsWT(wd)
            else C2 <- (1 - 2^(1 - beta))/(1 - 2^((1 - beta) * wd$
                  nlevels)) * sum(nsignal/p)
        }
        if(C2 < 0)
            stop("parameter C2 is negative")
        if(verbose == TRUE) cat("noise.level is: ", round(noise.level, 4), 
                "\nhyperparameters of the prior are: alpha = ", 
                alpha, "C1 = ", round(C1, 4), "beta = ", beta, 
                "C2 = ", round(C2, 4), "\n")    #   
#
# Bayesian Thresholding
#
        if(C1 == 0 | C2 == 0)
            wd <- nullevels(wd, levelstonu = levels)
        else {
            pr <- pmin(1, C2 * 2^( - beta * levels))
            rat <- tau2/(noise.level^2 + tau2)  #
            for(i in 1:nthresh) {
                d <- accessD(wd, level = levels[i])
                w <- (1 - pr[i])/pr[i]/sqrt((noise.level^2 * 
                  rat[i])/tau2[i]) * exp(( - rat[i] * d^2)/2/
                  noise.level^2)
                z <- 0.5 * (1 + pmin(w, 1))
                d <- sign(d) * pmax(0, rat[i] * abs(d) - 
                  noise.level * sqrt(rat[i]) * qnorm(z))
                wd <- putD(wd, level = levels[i], v = d)
            }
        }
        return(wd)
    }
    else if(policy == "cv") {
#
#
#       Threshold chosen by cross-validation
#
        if(verbose == TRUE) cat("Cross-validation policy\n")    #
        if(by.level == TRUE) stop(
                "Cross-validation policy does not permit by.level\n\t\t\tthresholding (yet)"
                )   #
#       Reconstruct the function for CWCV (this should be quick)
#
        ynoise <- wr(wd)
        thresh <- CWCV(ynoise = ynoise, x = 1:length(ynoise), 
            filter.number = wd$filter$filter.number, family = wd$
            filter$family, thresh.type = type, tol = cvtol, maxits=cvmaxits,
		verbose = 0, plot.it = FALSE, ll = min(levels))$xvthresh
        thresh <- rep(thresh, length = nthresh)
    }
    else if(policy == "fdr") {
#
#
#               Threshold chosen by FDR-procedure
#
        if(verbose == TRUE) cat("FDR policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh) {
                d <- c(d, accessD(wd, level = levels[i], 
                  boundary = boundary))
            }
            if(length(value) != 1)
                stop("Length of value should be 1")
            noise.level <- sqrt(dev(accessD(wd, level = (nlevelsWT(wd)-
                1))))
            minit <- length(d)
            dinit <- d
            thinit <- qnorm(1 - Q/2) * noise.level
            if(log(n, 2) > 12)
                ninit <- 3
            else {
                if(log(n, 2) > 10)
                  ninit <- 2
                else ninit <- 1
            }
            for(k in seq(1, ninit)) {
                dinit1 <- dinit[abs(dinit) >= thinit]
                minit <- length(dinit1)
                if(minit == 0)
                  thresh <- max(abs(d)) * 1.0001
                else {
                  thinit <- qnorm(1 - (Q * minit)/(2 * n)) * 
                    noise.level
                  minit1 <- length(dinit1[abs(dinit1) >= thinit
                    ])
                  if(minit1 == minit || minit1 == 0)
                    break
                  dinit <- dinit1
                }
            }
            if(noise.level > 0) {
                m <- length(d)
                minit <- length(dinit)
                p <- (2 - 2 * pnorm(abs(dinit)/noise.level))
                index <- order(p)
                j <- seq(1, minit)
                m0 <- max(j[p[index] <= (Q * j)/m])
                if(m0 != "NA" && m0 < minit)
                  thresh <- abs(dinit[index[m0]])
                else {
                  if(m0 == "NA")
                    thresh <- max(abs(dinit)) * 1.0001
                  else thresh <- 0
                }
            }
            else thresh <- 0
            thresh <- rep(thresh, length = nthresh)
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh[1], "\n", 
                  "sigma is: ", noise.level, "\n")
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for(i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = 
                  boundary)
                m <- length(d)
                noise.level <- sqrt(dev(d))
                thinit <- qnorm(1 - Q/2) * noise.level
                dinit <- d[abs(d) >= thinit]
                minit <- length(dinit)
                if(minit == 0)
                  thresh[i] <- max(abs(d)) * 1.0001
                else {
                  if(noise.level > 0) {
                    p <- (2 - 2 * pnorm(abs(dinit)/noise.level)
                      )
                    index <- order(p)
                    j <- seq(1, minit)
                    m0 <- max(j[p[index] <= (Q * j)/m])
                    if(m0 != "NA" && m0 < minit)
                      thresh[i] <- abs(dinit[index[m0]])
                    else {
                      if(m0 == "NA")
                        thresh[i] <- max(abs(dinit)) * 1.0001
                      else thresh[i] <- 0
                    }
                  }
                  else thresh[i] <- 0
                }
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], "is", 
                    thresh[i], "\n")
            }
        }
    }
    else if(policy == "op1") {
#
#
#       Ogden and Parzen's first policy
#
        if(verbose == TRUE) cat("Ogden and Parzen's first policy\n")
        if(by.level == FALSE)
            stop("Ogden and Parzen's first policy only computes level-dependent policies"
                )
        thresh <- TOthreshda1(ywd = wd, alpha = OP1alpha, verbose = 
            verbose, return.threshold = return.threshold)
        return(thresh)
    }
    else if(policy == "op2") {
#
#
#       Ogden and Parzen's second policy
#
        if(verbose == TRUE) cat("Ogden and Parzen's second policy\n")
        if(by.level == FALSE)
            stop("Ogden and Parzen's second policy only computes level-dependent policies"
                )
        thresh <- TOthreshda2(ywd = wd, alpha = OP1alpha, verbose = 
            verbose, return.threshold = return.threshold)
        return(thresh)
    }
    else if(policy == "manual") {
#
#
#       User supplied threshold policy
#
        if(verbose == TRUE) cat("Manual policy\n")
        thresh <- rep(value, length = nthresh)
        if(length(value) != 1 && length(value) != nthresh)
            warning("your threshold is not the same length as number of levels"
                )
    }
    else if(policy == "mannum") {
        if(verbose == TRUE) {
            cat("Manual policy using ", value, " of the")
            cat(" largest coefficients\n")
        }
        if(value < 1) {
            stop("Have to select an integer larger than 1 for value"
                )
        }
        else if(value > length(wd$D)) {
            stop(paste("There are only ", length(wd$D), 
                " coefficients, you specified ", value))
        }
        coefs <- wd$D
        scoefs <- sort(abs(coefs))
        scoefs <- min(rev(scoefs)[1:value])
        wd$D[abs(wd$D) < scoefs] <- 0
        return(wd)
    }
    else if(policy == "probability") {
#
#
#       Threshold is quantile based
#
        if(verbose == TRUE) cat("Probability policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh)
                d <- c(d, accessD(wd, level = levels[i], 
                  boundary = boundary))
            if(length(value) != 1)
                stop("Length of value should be 1")
            thresh <- rep(quantile(abs(d), prob = value), length = 
                nthresh)
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh[1], "\n")
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            if(length(value) == 1)
                value <- rep(value, nthresh)
            if(length(value) != nthresh)
                stop("Wrong number of probability values")
            for(i in 1:nthresh) {
                d <- accessD(wd, level = levels[i], boundary = 
                  boundary)
                thresh[i] <- quantile(abs(d), prob = value[i])
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
    }
    if(return.threshold == TRUE)
        return(thresh)
    for(i in 1:nthresh) {
        d <- accessD(wd, level = levels[i], boundary = boundary)
        if(type == "hard") {
            d[abs(d) <= thresh[i]] <- 0
        }
        else if(type == "soft") {
            d <- (d * (abs(d) - thresh[i]) * (abs(d) > thresh[i]))/
                abs(d)
            d[is.na(d)] <- 0
        }
        if(verbose == TRUE)
            cat("Level: ", levels[i], " there are ", sum(d == 0), 
                " zeroes\n")
        wd <- putD(wd, level = levels[i], v = d, boundary = boundary)
    }
    wd
}
"threshold.wd3D"<-
function(wd3D, levels = 3:(nlevelsWT(wd3D)- 1), type = "hard", policy = 
    "universal", by.level = FALSE, value = 0, dev = var, verbose = FALSE, 
    return.threshold = FALSE, ...)
{
    if(verbose == TRUE) cat("threshold.wd3D:\n")    #
#
#   Check class of wd3D
#
    if(verbose == TRUE)
        cat("Argument checking\n")
    ctmp <- class(wd3D)
    if(is.null(ctmp))
        stop("wd3D has no class")
    else if(ctmp != "wd3D")
        stop("wd3D is not of class wd3D")
    if(policy != "universal" && policy != "manual")
        stop("Only policys are universal, manual")
    if(type != "hard" && type != "soft")
        stop("Only hard or soft thresholding at  present")
    r <- range(levels)
    if(r[1] < 0)
        stop("levels out of range, level too small")
    if(r[2] > nlevelsWT(wd3D) - 1)
        stop(paste("levels out of range, level too big. Maximum level is ",
            nlevelsWT(wd3D) - 1))
    if(r[1] > nlevelsWT(wd3D) - 1) {
        warning("no thresholding done")
        return(wd3D)
    }
    if(r[2] < 0) {
        warning("no thresholding done")
        return(wd3D)
    }
    d <- NULL
    n <- (2^nlevelsWT(wd3D))^3
    nthresh <- length(levels)   #
#
#
#   Decide which policy to adopt
#       The next if-else construction should define a vector called
#       "thresh" that contains the threshold value for each level
#       in "levels". This may be the same threshold value
#       a global threshold.
#
    if(policy == "universal") {
#
#
#       Donoho and Johnstone's universal policy
#
        if(verbose == TRUE) cat("Universal policy...")
        if(by.level == FALSE) {
            if(verbose == TRUE)
                cat("All levels at once\n")
            for(i in 1:nthresh) {
                v <- accessD(wd3D, level = levels[i])
                d <- c(v$GHH, v$HGH, v$GGH, v$HHG, v$GHG, v$HGG,
                  v$GGG)
                if(levels[i] == 0)
                  d <- c(d, v$HHH)
            }
            noise.level <- sqrt(dev(d))
            nd <- length(d)
            thresh <- sqrt(2 * log(nd)) * noise.level
            if(verbose == TRUE)
                cat("Global threshold is: ", thresh, "\n")
            thresh <- rep(thresh, length = nthresh)
        }
        else {
            if(verbose == TRUE)
                cat("Level by level\n")
            thresh <- rep(0, length = nthresh)
            for(i in 1:nthresh) {
                v <- accessD(wd3D, level = levels[i])
                d <- c(v$GHH, v$HGH, v$GGH, v$HHG, v$GHG, v$HGG,
                  v$GGG)
                if(levels[i] == 0)
                  d <- c(d, v$HHH)
                noise.level <- sqrt(dev(d))
                nd <- length(d)
                thresh[i] <- sqrt(2 * log(nd)) * noise.level
                if(verbose == TRUE)
                  cat("Threshold for level: ", levels[i], 
                    " is ", thresh[i], "\n")
            }
        }
    }
    else if(policy == "manual") {
#
#
#       User supplied threshold policy
#
        if(verbose == TRUE) cat("Manual policy\n")
        thresh <- rep(value, length = nthresh)
        if(length(value) != 1 && length(value) != nthresh)
            warning("your threshold is not the same length as number of levels"
                )
    }
    if(return.threshold == TRUE)
        return(thresh)
    blocktypes <- c("GHH", "HGH", "GGH", "HHG", "GHG", "HGG", "GGG")
    for(i in 1:nthresh) {
        if(levels[i] == 0)
            lblocks <- c("HHH", blocktypes)
        else lblocks <- blocktypes
        nblocks <- length(lblocks)
        thedim <- rep(2^(levels[i]), 3)
        for(j in 1:nblocks) {
            d <- as.vector(accessD(wd3D, level = levels[i], block
                 = lblocks[j]))
            if(type == "hard") {
                d[abs(d) <= thresh[i]] <- 0
                if(verbose == TRUE)
                  cat("Level: ", levels[i], " there are ", sum(
                    d == 0), " zeroes\n")
            }
            else if(type == "soft") {
                d <- (d * (abs(d) - thresh[i]) * (abs(d) > 
                  thresh[i]))/abs(d)
                d[is.na(d)] <- 0
            }
            vinsert <- list(lev = levels[i], block = lblocks[j], a
                 = array(d, dim = thedim))
            wd3D <- putD(wd3D, v = vinsert)
        }
    }
    wd3D
}
"threshold.wp"<-
function(wp, levels = 3:(nlevelsWT(wp) - 1), dev = madmad, policy = "universal", 
    value = 0, by.level = FALSE, type = "soft", verbose = FALSE, return.threshold
     = FALSE, cvtol = 0.01, cvnorm = l2norm, add.history = TRUE, ...)
{
#
#   Do some arg checking
#
    rn <- range(levels)
    if(rn[1] < 0)
        stop("all selected levels must be larger than zero")
    if(rn[2] > nlevelsWT(wp) - 1)
        stop(paste("all selected levels must be smaller than", nlevelsWT(
            wp) - 1))
    nr <- nrow(wp$wp)
    nc <- ncol(wp$wp)   #
#
#   Figure out the threshold
#
    if(policy == "manual") {
        if(length(value) == 1) {
            if(verbose == TRUE)
                cat("Univariate threshold\n")
            threshv <- rep(value, length(levels))
        }
        else if(length(value) == length(levels)) {
            if(verbose == TRUE)
                cat("Multivariate threshold\n")
            threshv <- value
        }
        else stop("Manual policy. Your threshold vector is neither of length 1 or the length of the number of levels that you wish to threshold"
                )
    }
    else if(policy == "universal") {
        if(verbose == TRUE)
            cat("Universal threshold\n")
        if(by.level == FALSE) {
#
#       Global threshold
#
            d <- NULL
            for(lev in 1:length(levels)) {
                d <- c(d, accessD(wp, level = levels[lev]))
            }
            sigma <- dev(d)
            threshv <- sqrt(2 * log(nc) * sigma)
            threshv <- rep(threshv, length(levels))
        }
        else {
#
#
#       Level by level threshold
#
            threshv <- rep(0, length(levels))
            for(lev in 1:length(levels)) {
                d <- accessD(wp, level = levels[lev])
                sigma <- dev(d)
                threshv[lev] <- sqrt(2 * log(nc) * sigma)
            }
        }
    }
    if(verbose == TRUE) {
        cat("Threshold is ")
        print(threshv)
        cat("\n")
    }
#
#
#   Now apply the threshold
#
    if(return.threshold == TRUE)
        return(threshv)
    for(lev in 1:length(levels)) {
        if(verbose == TRUE) {
            cat("Applying threshold ", threshv[lev], " to level ", 
                levels[lev], "\n")
        }
        d <- accessD(wp, level = levels[lev])
        if(type == "hard")
            d[abs(d) <= threshv[lev]] <- 0
        else if(type == "soft")
            d <- sign(d) * (abs(d) - threshv[lev]) * (abs(d) > 
                threshv[lev])
        wp <- putD(wp, level = levels[lev], v = d)
    }
    wp$date <- c(wp$date, date())
    if(add.history == TRUE)
        wp$history <- c(wp$history, paste("Thresholded:", paste(
            as.character(threshv), collapse = "; "), "Levels: ", 
            paste(as.character(levels), collapse = "; "), 
            "Policy: ", policy, "Type: ", type))
    wp
}
"threshold.wst"<-
function(wst, levels = 3:(nlevelsWT(wst) - 1), dev = madmad, policy = "universal",
    value = 0, by.level = FALSE, type = "soft", verbose = FALSE, return.threshold
     = FALSE, cvtol = 0.01, cvnorm = l2norm, add.history = TRUE, ...)
{
#
#   Do some arg checking
#
    call <- match.call()
    rn <- range(levels)
    if(rn[1] < 0)
        stop("all selected levels must be larger than zero")
    if(rn[2] > nlevelsWT(wst) - 1)
        stop(paste("all selected levels must be smaller than", nlevelsWT(
            wst) - 1))
    nr <- nrow(wst$wp)
    nc <- ncol(wst$wp)  #
#
#   Figure out the threshold
#
    if(policy == "manual") {
        if(length(value) == 1) {
            if(verbose == TRUE)
                cat("Univariate threshold\n")
            threshv <- rep(value, length(levels))
        }
        else if(length(value) == length(levels)) {
            if(verbose == TRUE)
                cat("Multivariate threshold\n")
            threshv <- value
        }
        else stop("Manual policy. Your threshold vector is neither of length 1 or the length of the number of levels that you wish to threshold"
                )
    }
    else if(policy == "universal") {
        if(verbose == TRUE)
            cat("Universal threshold\n")
        if(by.level == FALSE) {
#
#       Global threshold
#
            d <- NULL
            for(lev in 1:length(levels)) {
                d <- c(d, accessD(wst, level = levels[lev]))
            }
            sigma <- dev(d)
            threshv <- sqrt(2 * log(nc) * sigma)
            threshv <- rep(threshv, length(levels))
        }
        else {
#
#
#       Level by level threshold
#
            threshv <- rep(0, length(levels))
            for(lev in 1:length(levels)) {
                d <- accessD(wst, level = levels[lev])
                sigma <- dev(d)
                threshv[lev] <- sqrt(2 * log(nc) * sigma)
            }
        }
    }
    else if(policy == "LSuniversal") {
        if(verbose == TRUE)
            cat("Local Spec universal threshold\n")
        if(by.level == FALSE) {
#
#       Global threshold
#
            d <- NULL
            for(lev in 1:length(levels)) {
                d <- c(d, accessD(wst, level = levels[lev]))
            }
            sigma <- dev(d)
            threshv <- log(nc) * sqrt(sigma)
            threshv <- rep(threshv, length(levels))
        }
        else {
#
#
#       Level by level threshold
#
            threshv <- rep(0, length(levels))
            for(lev in 1:length(levels)) {
                d <- accessD(wst, level = levels[lev])
                sigma <- dev(d)
                threshv[lev] <- log(nc) * sqrt(sigma)
            }
        }
    }
    else if(policy == "sure") {
        if(verbose == TRUE)
            cat("SURE threshold\n")
        if(by.level == FALSE) {
#
#       Global threshold
#
            d <- NULL
            for(lev in 1:length(levels)) {
                d <- c(d, accessD(wst, level = levels[lev]))
            }
            sigma <- sqrt(dev(d))
            threshv <- sigma * sure(d/sigma)
            threshv <- rep(threshv, length(levels))
        }
        else {
#
#
#       Level by level threshold
#
            threshv <- rep(0, length(levels))
            for(lev in 1:length(levels)) {
                d <- accessD(wst, level = levels[lev])
                sigma <- sqrt(dev(d))
                threshv[lev] <- sigma * sure(d/sigma)
            }
        }
    }
    else if(policy == "cv") {
        if(verbose == TRUE)
            cat("Cross-validation threshold\n")
        ynoise <- AvBasis(wst)
        if(by.level == TRUE) {
            if(verbose == TRUE)
                cat("by-level\n")
            if(length(levels) != 1)
                warning(
                  "Taking minimum level as first level for level-dependent cross-validation"
                  )
            levels <- min(levels):(nlevelsWT(wst) - 1)
            threshv <- wstCVl(ndata = ynoise, ll = min(levels), 
                type = type, filter.number = wst$filter$
                filter.number, family = wst$filter$family, tol
                 = cvtol, verbose = 0, plot.it = FALSE, norm = 
                cvnorm, InverseType = "average")$xvthresh
            if(verbose == TRUE)
                cat("Cross-validation threshold is ", threshv, 
                  "\n")
        }
        else {
            if(verbose == TRUE)
                cat("global\n")
            threshv <- wstCV(ndata = ynoise, ll = min(levels), type
                 = type, filter.number = wst$filter$
                filter.number, family = wst$filter$family, tol
                 = cvtol, verbose = 0, plot.it = FALSE, norm = 
                cvnorm, InverseType = "average")$xvthresh
            threshv <- rep(threshv, length(levels))
        }
    }
    else {
        stop(paste("Unknown policy: ", policy))
    }
    if(verbose == TRUE) {
        cat("Threshold is ")
        print(threshv)
        cat("\n")
    }
#
#
#   Now apply the threshold
#
    if(return.threshold == TRUE)
        return(threshv)
    for(lev in 1:length(levels)) {
        if(verbose == TRUE) {
            cat("Applying threshold ", threshv[lev], " to level ", 
                levels[lev], "(type is ", type, ")\n")
        }
        d <- accessD(wst, level = levels[lev])
        if(type == "hard")
            d[abs(d) <= threshv[lev]] <- 0
        else if(type == "soft")
            d <- sign(d) * (abs(d) - threshv[lev]) * (abs(d) > 
                threshv[lev])
        wst <- putD(wst, level = levels[lev], v = d)
    }
    wst$date <- c(wst$date, date())
    if(add.history == TRUE)
        wst$history <- c(wst$history, paste("Thresholded:", paste(
            as.character(threshv), collapse = "; "), "Levels: ", 
            paste(as.character(levels), collapse = "; "), 
            "Policy: ", policy, "Type: ", type))
    wst
}
"tpwd"<-
function(image, filter.number = 10, family = "DaubLeAsymm", verbose = FALSE)
{
    if(!is.matrix(image))
        stop("image should be a matrix")
    nr <- nrow(image)
    lr <- IsPowerOfTwo(nr)
    if(is.na(lr))
        stop(paste("Number of rows (", nr, ") should be a power of 2.")
            )
    nc <- ncol(image)
    lc <- IsPowerOfTwo(nc)
    if(is.na(lc))
        stop(paste("Number of cols (", nc, ") should be a power of 2.")
            )
    bc <- "periodic"
    type <- "wavelet"
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2)
    if(is.null(nbc))
        stop("Unknown boundary condition")
    ntype <- switch(type,
        wavelet = 1,
        station = 2)    #
#
# Select the appropriate filter
#
    if(verbose == TRUE)
        cat("...done\nFilter...")
    filter <- filter.select(filter.number = filter.number, family = family)
        #
#
# Build the first/last database
#
    if(verbose == TRUE)
        cat("...selected\nFirst/last database...")
    fl.dbaseR <- first.last(LengthH = length(filter$H), DataLength = nr, 
        type = type, bc = bc)   #
    fl.dbaseC <- first.last(LengthH = length(filter$H), DataLength = nc, 
        type = type, bc = bc)   #
    error <- 0
    answer <- .C("tpwd",
        image = as.double(image),
        nr = as.integer(nr),
        nc = as.integer(nc),
        lr = as.integer(lr),
        lc = as.integer(lc),
        firstCr = as.integer(fl.dbaseR$first.last.c[, 1]),
        lastCr = as.integer(fl.dbaseR$first.last.c[, 2]),
        offsetCr = as.integer(fl.dbaseR$first.last.c[, 3]),
        firstDr = as.integer(fl.dbaseR$first.last.d[, 1]),
        lastDr = as.integer(fl.dbaseR$first.last.d[, 2]),
        offsetDr = as.integer(fl.dbaseR$first.last.d[, 3]),
        firstCc = as.integer(fl.dbaseC$first.last.c[, 1]),
        lastCc = as.integer(fl.dbaseC$first.last.c[, 2]),
        offsetCc = as.integer(fl.dbaseC$first.last.c[, 3]),
        firstDc = as.integer(fl.dbaseC$first.last.d[, 1]),
        lastDc = as.integer(fl.dbaseC$first.last.d[, 2]),
        offsetDc = as.integer(fl.dbaseC$first.last.d[, 3]),
        ntype = as.integer(ntype),
        nbc = as.integer(nbc),
        H = as.double(filter$H),
        LengthH = as.integer(length(filter$H)),
        error = as.integer(error), PACKAGE = "wavethresh")
    theanswer <- list(tpwd = matrix(answer$image, nrow = nr, ncol = nc), 
        filter.number = filter.number, family = family, type = type, bc
         = bc, date = date())
    class(theanswer) <- "tpwd"
    theanswer
}
"tpwr"<-
function(tpwdobj, verbose = FALSE)
{
    if(class(tpwdobj) != "tpwd")
        stop("tpwdobj is not of class tpwd")
    nr <- nrow(tpwdobj$tpwd)
    lr <- IsPowerOfTwo(nr)
    nc <- ncol(tpwdobj$tpwd)
    lc <- IsPowerOfTwo(nc)
    bc <- tpwdobj$bc
    type <- tpwdobj$type
    nbc <- switch(bc,
        periodic = 1,
        symmetric = 2)
    ntype <- switch(type,
        wavelet = 1,
        station = 2)    #
#
# Select the appropriate filter
#
    if(verbose == TRUE)
        cat("...done\nFilter...")
    filter <- filter.select(filter.number = tpwdobj$filter.number, family
         = tpwdobj$family)  #
#
# Build the first/last database
#
    if(verbose == TRUE)
        cat("...selected\nFirst/last database...")
    fl.dbaseR <- first.last(LengthH = length(filter$H), DataLength = nr, 
        type = type, bc = bc)   #
    fl.dbaseC <- first.last(LengthH = length(filter$H), DataLength = nc, 
        type = type, bc = bc)   #
    error <- 0
    answer <- .C("tpwr",
        image = as.double(tpwdobj$tpwd),
        nr = as.integer(nr),
        nc = as.integer(nc),
        lr = as.integer(lr),
        lc = as.integer(lc),
        firstCr = as.integer(fl.dbaseR$first.last.c[, 1]),
        lastCr = as.integer(fl.dbaseR$first.last.c[, 2]),
        offsetCr = as.integer(fl.dbaseR$first.last.c[, 3]),
        firstDr = as.integer(fl.dbaseR$first.last.d[, 1]),
        lastDr = as.integer(fl.dbaseR$first.last.d[, 2]),
        offsetDr = as.integer(fl.dbaseR$first.last.d[, 3]),
        firstCc = as.integer(fl.dbaseC$first.last.c[, 1]),
        lastCc = as.integer(fl.dbaseC$first.last.c[, 2]),
        offsetCc = as.integer(fl.dbaseC$first.last.c[, 3]),
        firstDc = as.integer(fl.dbaseC$first.last.d[, 1]),
        lastDc = as.integer(fl.dbaseC$first.last.d[, 2]),
        offsetDc = as.integer(fl.dbaseC$first.last.d[, 3]),
        ntype = as.integer(ntype),
        nbc = as.integer(nbc),
        H = as.double(filter$H),
        LengthH = as.integer(length(filter$H)),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(answer$error != 0)
        stop(paste("Error code was ", answer$error))
    theanswer <- matrix(answer$image, nrow = nr, ncol = nc)
    theanswer
}
"uncompress"<-
function(...)
UseMethod("uncompress")
"uncompress.default"<-
function(v, verbose = FALSE, ...)
{
    ctmp <- class(v)
    if(is.null(ctmp)) {
        stop("Object v has no class")
    }
    else if(ctmp == "uncompressed") {
        if(verbose == TRUE)
            cat("Not compressed\n")
        return(unclass(v$vector))
    }
    else if(ctmp == "compressed") {
        answer <- rep(0, length = v$original.length)
        answer[v$position] <- v$values
        if(verbose == TRUE)
            cat("Uncompressed to length ", length(answer), "\n")
        return(answer)
    }
    else stop("v has unknown class")
}
"uncompress.imwdc"<-
function(x, verbose = FALSE, ...)
{
    if(verbose == TRUE)
        cat("Argument checking\n")
    ctmp <- class(x)
    if(is.null(ctmp))
        stop("imwd has no class")
    else if(ctmp != c("imwdc"))
        stop("imwd is not of class imwdc")
    unsquished <- list(nlevels = nlevelsWT(x), fl.dbase = x$fl.dbase, 
        filter = x$filter, w0Lconstant = x$w0Lconstant, bc = x$
        bc, type = x$type)   #
#
#       Go round loop compressing each set of coefficients
#
    for(level in 0:(nlevelsWT(x)- 1)) {
        if(verbose == TRUE)
            cat("Level ", level, "\n\t")
        nm <- lt.to.name(level, "CD")
        if(verbose == TRUE)
            cat("CD\t")
        unsquished[[nm]] <- uncompress.default(x[[nm]], verbose = 
            verbose)
        nm <- lt.to.name(level, "DC")
        if(verbose == TRUE)
            cat("\tDC\t")
        unsquished[[nm]] <- uncompress.default(x[[nm]], verbose = 
            verbose)
        nm <- lt.to.name(level, "DD")
        if(verbose == TRUE)
            cat("\tDD\t")
        unsquished[[nm]] <- uncompress.default(x[[nm]], verbose = 
            verbose)
    }
    class(unsquished) <- "imwd"
    if(verbose == TRUE)
        cat("Overall inflation: Was: ", w <- object.size(x), " Now:",
            s <- object.size(unsquished), " (", signif((100 * s)/w, 
            digits=3), "%)\n")
    unsquished
}
"wavegrow"<-
function(n = 64, filter.number = 10, family = "DaubLeAsymm", type = "wavelet", 
    random = TRUE, read.value = TRUE, restart = FALSE)
{
    nlev <- IsPowerOfTwo(n)
    if(is.na(nlev))
        stop("n is not a power of two")
    coords <- vector("list", nlev)
    if(type == "wavelet") {
        x <- 1:(n/2)
        coords[[nlev]] <- x
        nn <- n/2
        for(i in (nlev - 1):1) {
            x1 <- x[seq(1, nn - 1, 2)]
            x2 <- x[seq(2, nn, 2)]
            x <- (x1 + x2)/2
            nn <- nn/2
            coords[[i]] <- x
        }
    }
    else for(i in 1:nlev)
            coords[[i]] <- 1:n
    if(is.null(dev.list()))
        stop("Please start 2 graphical devices before using me")
    if(length(dev.list()) < 2)
        stop("Please start another graphics device\n")
    ndev <- length(dev.list())
    gd1 <- dev.list()[ndev - 1]
    gd2 <- dev.list()[ndev]
    v <- rnorm(n, sd = 1e-10)
    vwr <- v
    vwdS <- wd(v, filter.number = filter.number, family = family, type = type)
    toplev <- nlevelsWT(vwdS) - 1
    ans <- "y"
    while(ans == "y" | ans == "yes" | ans == "Y") {
        dev.set(which = gd1)
        ts.plot(v)
        dev.set(which = gd2)
        plot(vwdS, NotPlotVal = 0)
        while(1) {
            co <- locator(1)
            if(is.null(co))
                break
            lev <- 1 + toplev - round(co$y)
            cvec <- coords[[lev + 1]]
            ix <- (cvec - co$x)^2
            nvec <- length(cvec)
            ix <- (1:nvec)[ix == min(ix)]
            if(type == "station") {
                ix <- ix - 2^(nlev - lev - 1)
                ix <- ((ix - 1) %% n) + 1
            }
            cat("Level ", lev, " Coordinate ", ix, "\n")
            if(random == TRUE)
                new <- rnorm(1)
            else {
                if(read.value == TRUE) {
                  cat("Type in coefficient value ")
                  new <- scan(n = 1)
                }
                else new <- 1
            }
            v <- accessD(vwdS, lev = lev)
            v[ix] <- new
            vwdS <- putD(vwdS, lev = lev, v = v)
            plot(vwdS, NotPlotVal = 0)
            dev.set(which = gd1)
            if(type == "station") {
                vwdWST <- convert(vwdS)
                vwr <- AvBasis(vwdWST)
            }
            else vwr <- wr(vwdS)
            ts.plot(vwr)
            dev.set(which = gd2)
            if(restart == TRUE) {
                v <- rep(1, n)
                vwdS <- wd(v, filter.number = filter.number, family = 
                  family, type = type)
            }
        }
        cat("Do you want to continue? ")
        ans <- readline()
        if(ans == "y" | ans == "yes" | ans == "Y") {
            v <- rnorm(n, sd = 1e-10)
            vwdS <- wd(v, filter.number = filter.number, family = family, 
                type = type)
        }
    }
    return(list(ts = vwr, wd = vwdS))
}
"wd.int"<-
function(data, preferred.filter.number, min.scale, precond)
{
    storage.mode(data) <- "double"
    storage.mode(preferred.filter.number) <- "integer"
    storage.mode(min.scale) <- "integer"
    storage.mode(precond) <- "logical"
    size <- length(data)
    storage.mode(size) <- "integer"
    max.scale <- log(size, 2)
    filter.history <- integer(max.scale - min.scale)
    temp <- .C("dec",
        vect = data,
        size,
        preferred.filter.number,
        min.scale,
        precond,
        history = filter.history, PACKAGE = "wavethresh")
    wav.int.object <- list(transformed.vector = temp$vect, current.scale = 
        min.scale, filters.used = temp$history, preconditioned = 
        precond, date = date())
    return(wav.int.object)
}
"wd3D"<-
function(a, filter.number = 10, family = "DaubLeAsymm")
{
    d <- dim(a)
    if(length(d) != 3)
        stop(paste("a is not a three-dimensional object"))
    for(i in 1:3)
        if(is.na(IsPowerOfTwo(d[i]))) stop(paste("Dimension ", i, 
                " of a is not of dyadic length"))
    if(any(d != d[1]))
        stop("Number of elements in each dimension is not identical")
    error <- 0
    nlevels <- IsPowerOfTwo(d[1])
    H <- filter.select(filter.number = filter.number, family = family)$H
    ans <- .C("wd3D",
        Carray = as.double(a),
        size = as.integer(d[1]),
        H = as.double(H),
        LengthH = as.integer(length(H)),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(ans$error != 0)
        stop(paste("Error code was ", ans$error))
    l <- list(a = array(ans$Carray, dim = d), filter.number = filter.number,
        family = family, date = date(), nlevels = nlevels)
    class(l) <- "wd3D"
    l
}
"wp"<-
function(data, filter.number = 10, family = "DaubLeAsymm", verbose = FALSE)
{
    if(verbose == TRUE)
        cat("Argument checking...")
    DataLength <- length(data)  #
#
# Check that we have a power of 2 data elements
#
    nlevels <- log(DataLength)/log(2)
    if(round(nlevels) != nlevels)
        stop("The length of data is not a power of 2")  #
    if(verbose == TRUE) {
        cat("There are ", nlevels, " levels\n")
    }
#
# Select the appropriate filter
#
    if(verbose == TRUE)
        cat("...done\nFilter...")
    filter <- filter.select(filter.number = filter.number, family = family)
        #
#
# Compute the decomposition
#
    if(verbose == TRUE)
        cat("Decomposing...\n")
    newdata <- c(rep(0, DataLength * nlevels), data)
    wavelet.packet <- .C("wavepackde",
        newdata = as.double(newdata),
        DataLength = as.integer(DataLength),
        levels = as.integer(nlevels),
        H = as.double(filter$H),
        LengthH = as.integer(length(filter$H)), PACKAGE = "wavethresh")
    wpm <- matrix(wavelet.packet$newdata, ncol = DataLength, byrow = TRUE)
    wp <- list(wp = wpm, nlevels = nlevels, filter = filter, date = date())
    class(wp) <- "wp"
    wp
}
"wpst"<-
function(data, filter.number = 10, family = "DaubLeAsymm", FinishLevel = 0)
{
    nlev <- nlevelsWT(data)
    n <- length(data)
    if(FinishLevel < 0)
        stop("FinishLevel must be larger than zero")
    else if(FinishLevel >= nlev)
        stop(paste("FinishLevel must be < ", nlev)) #   
    lansvec <- n * (2 * n - 1)
    ansvec <- rep(0, lansvec)   #
#
#   Now create vector that keeps track of where levels start/stop
#
#   Note that the vector avixstart stores index entry values in C
#   notation. If you use it in Splus you'll have to add on 1
#
    npkts <- function(level, nlev)
    4^(nlev - level)
    pktlength <- function(level)
    2^level
    avixstart <- rep(0, nlev + 1)
    for(i in 1:nlev)
        avixstart[i + 1] <- avixstart[i] + npkts(i - 1, nlev) * 
            pktlength(i - 1)    #
#
#   Copy in original data
#
    ansvec[(avixstart[nlev + 1] + 1):lansvec] <- data   #
#
#   Call the C routine
#
    filter <- filter.select(filter.number = filter.number, family = family)
    ans <- .C("wpst",
        ansvec = as.double(ansvec),
        lansvec = as.integer(lansvec),
        nlev = as.integer(nlev),
        FinishLevel = as.integer(FinishLevel),
        avixstart = as.integer(avixstart),
        H = as.double(filter$H),
        LengthH = as.integer(length(filter$H)),
        error = as.integer(0), PACKAGE = "wavethresh")
    rv <- list(wpst = ans$ansvec, nlevels = nlev, avixstart = avixstart, 
        filter = filter, date = date())
    class(rv) <- "wpst"
    rv
}
"wpst2discr"<-
function(wpstobj, groups)
{
#
#   Function to convert wpst object and associated groups vector into
#   data matrix and k vector required as the input to the discr function.
#
#   Input:  wpstobj: a wpst object of a time-series
#       groups: a vector of length ncases containing the group
#           membership of each case.
#
#   Returns: wpstm  - a matrix. Number of rows is the number of cases
#           The rows are ordered according to the group
#           memberships of the cases. E.g. The first n1 rows
#           contain the group 1 cases, the second n2 rows
#           contain the group 2 cases, ... the ng rows
#           contain the group g cases.
#
#       level   - a vector of length npkts. Each entry refers to
#           the level that the col of wpstm comes from.
#
#       pktix   - a vector of length npkts. Each entry refers to
#           the packet index that the col of wpstm comes from.
#
#
#       k   - a vector of length ng (the number of groups).
#            k[1] contains the number of members for group 1, 
#            k[2] contains the number of members for group 2, ...
#            k[ng] contains the number of members for group ng.
#
#
#
    J <- nlev <- nlevelsWT(wpstobj)
    grot <- compgrot(J, filter.number=2)
    nbasis <- 2 * (2^nlev - 1)
    ndata <- 2^nlev
    m <- matrix(0, nrow = ndata, ncol = nbasis)
    level <- rep(0, nbasis)
    pktix <- rep(0, nbasis)
    cnt <- 1
    cat("Level: ")
    for(j in 0:(nlev - 1)) {
        cat(j, " ")
        lcnt <- 0
        npkts <- 2^(nlev - j)
        prcnt <- as.integer(npkts/10)
	if (prcnt == 0)
		prcnt <- 1
        for(i in 0:(npkts - 1)) {
            pkcoef <- guyrot(accessD(wpstobj, level = j, index = i),
                grot[J - j])/(sqrt(2)^(J - j))
            m[, cnt] <- log(pkcoef^2)
            level[cnt] <- j
            pktix[cnt] <- i
            lcnt <- lcnt + 1
            cnt <- cnt + 1
            if(lcnt %% prcnt == 0) {
                lcnt <- 0
                cat(".")
            }
        }
        cat("\n")
    }
    cat("\n")
    l <- list(m = m, groups = groups, level = level, pktix = pktix, nlevels = J)
    class(l) <- "w2d"
    l
}
"wpstCLASS"<-
function(newTS, wpstDO)
{
#
#
# Apply wpst to new TS
#
    newwpst <- wpst(newTS, filter.number = wpstDO$filter$filter.number, 
        family = wpstDO$filter$family)  #
#
# Extract the "best packets"
#
    goodlevel <- wpstDO$BP$level
    goodpkt <- wpstDO$BP$pkt
    npkts <- length(goodpkt)
    ndata <- length(newTS)
    m <- matrix(0, nrow = ndata, ncol = npkts)
    J <- nlevelsWT(newwpst)
    grot <- compgrot(J, filter.number=2)
    for(i in 1:npkts) {
        j <- goodlevel[i]
        m[, i] <- guyrot(accessD(newwpst, level = j, index = goodpkt[i]
            ), grot[J - j])/(sqrt(2)^(J - j))
        m[, i] <- log(m[, i]^2)
    }
    mTd <- predict(wpstDO$BPd$dm, m)

	l <- list(BasisMatrix=m, BasisMatrixDM=m%*%wpstDO$BPd$dm$scaling,
		wpstDO=wpstDO, PredictedOP=mTd, PredictedGroups=mTd$class)
	class(l) <- "wpstCL"
	l
}
"wr"<-
function(...)
UseMethod("wr")
"wr.int"<-
function(wav.int.object, ...)
{
    data <- wav.int.object$transformed.vector
    storage.mode(data) <- "double"
    size <- length(data)
    storage.mode(size) <- "integer"
    filter.history <- wav.int.object$filters.used
    storage.mode(filter.history) <- "integer"
    current.scale <- wav.int.object$current.scale
    storage.mode(current.scale) <- "integer"
    precond <- wav.int.object$preconditioned
    storage.mode(precond) <- "logical"
    temp <- .C("rec",
        vect = data,
        size,
        filter.history,
        current.scale,
        precond, PACKAGE = "wavethresh")
    return(temp$vect)
}
"wr.mwd"<-
function(...)
{
#calling mwr directly would be better but
#just in case...
    mwr(...)
}
"wr3D"<-
function(obj)
{
    ClassObj <- class(obj)
    if(is.null(ClassObj))
        stop("obj has no class")
    if(ClassObj != "wd3D")
        stop("obj is not of class wd3D")
    Carray <- obj$a
    H <- filter.select(filter.number = obj$filter.number, family = obj$
        family)$H
    answer <- .C("wr3D",
        Carray = as.double(Carray),
        truesize = as.integer(dim(Carray)[1]),
        H = as.double(H),
        LengthH = as.integer(length(H)),
        error = as.integer(0), PACKAGE = "wavethresh")
    array(answer$Carray, dim = dim(Carray))
}
"wst2D"<-
function(m, filter.number = 10, family = "DaubLeAsymm")
{
    nr <- nrow(m)
    J <- IsPowerOfTwo(nr)
    dimv <- c(J, 2 * nr, 2 * nr)
    am <- array(0, dim = dimv)
    filter <- filter.select(filter.number = filter.number, family = family)
    error <- 0
    ans <- .C("SWT2Dall",
        m = as.double(m),
        nm = as.integer(nr),
        am = as.double(am),
        J = as.integer(J),
        H = as.double(filter$H),
        LengthH = as.integer(length(filter$H)),
        error = as.integer(error), PACKAGE = "wavethresh")
    if(ans$error != 0)
        stop(paste("Error code was ", ans$error))
    l <- list(wst2D = array(ans$am, dim = dimv), nlevels = J, filter = 
        filter, date = date())
    class(l) <- "wst2D"
    l
}
"wstCV"<-
function(ndata, ll = 3, type = "soft", filter.number = 10, family = 
    "DaubLeAsymm", tol = 0.01, verbose = 0, plot.it = FALSE, norm = l2norm, 
    InverseType = "average", uvdev = madmad)
{
    nlev <- log(length(ndata))/log(2)
    levels <- ll:(nlev - 1)
    nwst <- wst(ndata, filter.number = filter.number, family = family)
    uv <- threshold(nwst, levels = levels, type = type, policy = 
        "universal", dev = madmad, return.thresh = TRUE)[1]
    if(verbose == 1)
        cat("Now optimising cross-validated error estimate\n")
    levels <- ll:(nlev - 2)
    R <- 0.61803399000000003
    C <- 1 - R
    ax <- 0
    bx <- uv/2
    cx <- uv
    x0 <- ax
    x3 <- cx
    if(abs(cx - bx) > abs(bx - ax)) {
        x1 <- bx
        x2 <- bx + C * (cx - bx)
    }
    else {
        x2 <- bx
        x1 <- bx - C * (bx - ax)
    }
    fa <- GetRSSWST(ndata, threshold = ax, levels = levels, type = type, 
        filter.number = filter.number, family = family, norm = norm, 
        verbose = verbose, InverseType = InverseType)
    cat("Done 1\n")
    fb <- GetRSSWST(ndata, threshold = bx, levels = levels, type = type, 
        filter.number = filter.number, family = family, norm = norm, 
        verbose = verbose, InverseType = InverseType)
    cat("Done 2\n")
    fc <- GetRSSWST(ndata, threshold = cx, levels = levels, type = type, 
        filter.number = filter.number, family = family, norm = norm, 
        verbose = verbose, InverseType = InverseType)
    cat("Done 3\n")
    f1 <- GetRSSWST(ndata, threshold = x1, levels = levels, type = type, 
        filter.number = filter.number, family = family, norm = norm, 
        verbose = verbose, InverseType = InverseType)
    cat("Done 4\n")
    f2 <- GetRSSWST(ndata, threshold = x2, levels = levels, type = type, 
        filter.number = filter.number, family = family, norm = norm, 
        verbose = verbose, InverseType = InverseType)
    cat("Done 5\n")
    xkeep <- c(ax, cx, x1, x2)
    fkeep <- c(fa, fc, f1, f2)
    if(plot.it == TRUE) {
        plot(c(ax, bx, cx), c(fa, fb, fc))
        text(c(x1, x2), c(f1, f2), lab = c("1", "2"))
    }
    cnt <- 3
    while(abs(x3 - x0) > tol * (abs(x1) + abs(x2))) {
        if(verbose > 0) {
            cat("x0=", x0, "x1=", x1, "x2=", x2, "x3=", x3, "\n")
            cat("f1=", f1, "f2=", f2, "\n")
        }
        if(f2 < f1) {
            x0 <- x1
            x1 <- x2
            x2 <- R * x1 + C * x3
            f1 <- f2
            f2 <- GetRSSWST(ndata, threshold = x2, levels = levels, 
                type = type, filter.number = filter.number, 
                family = family, norm = norm, verbose = verbose,
                InverseType = InverseType)
            if(verbose == 2) {
                cat("SSQ: ", signif(f2, digits=3), "\n")
            }
            else if(verbose == 1)
                cat(".")
            xkeep <- c(xkeep, x2)
            fkeep <- c(fkeep, f2)
            if(plot.it == TRUE)
                text(x2, f2, lab = as.character(cnt))
            cnt <- cnt + 1
        }
        else {
            x3 <- x2
            x2 <- x1
            x1 <- R * x2 + C * x0
            f2 <- f1
            f1 <- GetRSSWST(ndata, threshold = x1, levels = levels, 
                type = type, filter.number = filter.number, 
                family = family, norm = norm, verbose = verbose,
                InverseType = InverseType)
            if(verbose == 2)
                cat("SSQ: ", signif(f1, digits=3), "\n")
            else if(verbose == 1)
                cat(".")
            xkeep <- c(xkeep, x1)
            fkeep <- c(fkeep, f1)
            if(plot.it == TRUE)
                text(x1, f1, lab = as.character(cnt))
            cnt <- cnt + 1
        }
    }
    if(f1 < f2)
        tmp <- x1
    else tmp <- x2
    x1 <- tmp/sqrt(1 - log(2)/log(length(ndata)))
    if(verbose == 1)
        cat("Correcting to ", x1, "\n")
    else if(verbose == 1)
        cat("\n")
    g <- sort.list(xkeep)
    xkeep <- xkeep[g]
    fkeep <- fkeep[g]
    if(verbose >= 1) {
        cat("Reconstructing CV \n")
    }
    nwstT <- threshold(nwst, type = type, levels = levels, policy = 
        "manual", value = x1)   #
#
#   Now threshold the top level using universal thresholding
#
    nwstT <- threshold(nwstT, type = type, levels = nlevelsWT(nwstT) - 1, 
        policy = "universal", dev = uvdev)
    xvwr <- AvBasis.wst(nwstT)
    list(ndata = ndata, xvwr = xvwr, xvwrWSTt = nwstT, uvt = uv, xvthresh
         = x1, xkeep = xkeep, fkeep = fkeep)
}
"wstCVl"<-
function(ndata, ll = 3, type = "soft", filter.number = 10, family = 
    "DaubLeAsymm", tol = 0.01, verbose = 0, plot.it = FALSE, norm = l2norm, 
    InverseType = "average", uvdev = madmad)
{
    nlev <- log(length(ndata))/log(2)
    levels <- ll:(nlev - 2)
    nwst <- wst(ndata, filter.number = filter.number, family = family)
    uv <- threshold(nwst, levels = levels, type = type, policy = 
        "universal", dev = madmad, return.thresh = TRUE)[1]
    if(verbose == 1)
        cat("Now optimising cross-validated error estimate\n")
    upper <- rep(uv, length(levels))
    lower <- rep(0, length(levels))
    start <- (lower + upper)/2
    answer <- nlminb(start = start, objective = wvcvlrss, lower = lower, 
        upper = upper, ndata = ndata, levels = levels, type = type, 
        filter.number = filter.number, family = family, norm = norm, 
        verbose = verbose, InverseType = InverseType, control = list(rel.tol = tol))
    x1 <- answer$par
    if(verbose >= 2)
        thverb <- TRUE
    else thverb <- FALSE
    xvwrWSTt <- threshold.wst(nwst, levels = levels, policy = "manual", 
        value = x1, verbose = thverb)   #
#       Now threshold the top level using universal thresholding
#
    lastuvt <- threshold(xvwrWSTt, type = type, levels = nlevelsWT(xvwrWSTt) - 
        1, policy = "universal", dev = uvdev, return.thresh = TRUE)
    xvwrWSTt <- threshold(xvwrWSTt, type = type, levels = nlevelsWT(xvwrWSTt) -
        1, policy = "manual", value = lastuvt)
    xvwr <- AvBasis.wst(xvwrWSTt)
    list(ndata = ndata, xvwr = xvwr, xvwrWSTt = xvwrWSTt, uvt = uv, 
        xvthresh = c(x1, lastuvt), optres = answer)
}
"wvcvlrss"<-
function(threshold, ndata, levels, type, filter.number, family, norm, verbose, 
    InverseType)
{
    answer <- GetRSSWST(ndata = ndata, threshold = threshold, levels = 
        levels, family = family, filter.number = filter.number, type = 
        type, norm = norm, verbose = verbose, InverseType = InverseType
        )
    return(answer)
}
"wvmoments"<-
function(filter.number = 10, family = "DaubLeAsymm", moment = 0, 
    scaling.function = FALSE)
{
    WV <- draw.default(filter.number = filter.number, family = family, 
        plot.it = FALSE, enhance = FALSE, resolution = 32768, scaling.function = 
        scaling.function)
    intfn <- function(x, moment, xwv, ywv)
    {
        x^moment * approx(x = xwv, y = ywv, xout = x, rule = 2)$y
    }
    plot(WV$x, intfn(WV$x, moment = moment, WV$x, WV$y), type = "l")
    integrate(intfn, lower = -7, upper = 7, moment = moment, xwv = WV$x, 
        ywv = WV$y, subdivisions = 1000, keep.xy = TRUE)
}
"wvrelease"<-
function()
{
    packageStartupMessage("WaveThresh: R wavelet software, release 4.6.6, installed\n")
    packageStartupMessage("Copyright Guy Nason and others 1993-2013\n")
    packageStartupMessage("Note: nlevels has been renamed to nlevelsWT\n")
}
