l1median <- function (X, MaxStep = 200, ItTol = 10^-8, trace = 0, m.init = .colMedians (X))
{
	if (is.null (dim(X)))
		return (median (X))
	l1median_NLM (X = X, maxit = MaxStep, tol = ItTol, trace = trace, m.init = m.init)$par
}

#l1median = function (X, MaxStep = 200, ItTol = 10^-8, trace = 0)
#{
#	if (trace >= 0)
#		warning ("This function (pcaPP::l1median)is outdated.\r\nFor better performance try any of pcaPP::l1median_* instead. Preferably pcaPP::l1median_NLM.\r\nOtherwise use (trace = -1) for suppressing this warning. ")
#
#	if (class (X) != "matrix")
#	{
#		if (class (X) == "data.frame")
#			X = as.matrix(X)
#		else
#			X = matrix(X, ncol = 1)
#	}
#
#	ret = .C ("l1median", PACKAGE="pcaPP",
#		as.double (X),
#		as.integer (nrow(X)),
#		as.integer (ncol(X)),
#		med = double (ncol(X)),
#		ret = integer(1),
#		as.integer (MaxStep),
#		as.double (ItTol)
#		)
#
#	if (ret$ret != 0)
#		return (ret$med)
#	stop("iteration failed")
#}


l1median_BFGS <-
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = .colMedians (X), REPORT = 10, ...)
{
	X <- .Conv2Matrix (X, substitute (X))

	if (length (m.init) != ncol (X))
		stop (paste ("length of vector m.init (=", length (m.init), ") does not match the number of columns of data object X (=", ncol (X),")", sep = ""))
		
		
	ret = .C ("l1median_BFGS", PACKAGE="pcaPP", NAOK = TRUE, 
		par = as.integer (c(dim (X), maxit, trace, REPORT)),
		npar.out = integer (4),
		dpar = as.double (c(-Inf, tol)),
		dpar.out = double (1),
		as.double (X),
		#as.double (pscale),
		med = as.double (m.init)#double (ncol(X))
		)

	if (trace >= 1)
		cat ("l1median returned", ret$npa, ret$dpar, "\r\n") ;

	return (list (par = ret$med, value = ret$dpar.out[1], code = ret$npar.out [1], iterations = ret$npar.out [2], iterations_gr = ret$npar.out [3], time = ret$npar.out[4]))
}


l1median_CG <-
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = .colMedians (X), ...)#, type = 1)
{
	type = 1
	
	X <- .Conv2Matrix (X, substitute (X))

	if (length (m.init) != ncol (X))
		stop (paste ("length of vector m.init (=", length (m.init), ") does not match the number of columns of data object X (=", ncol (X),")", sep = ""))
		
	if (type < 1 || type > 3)
		stop ("parameter type MUST be either 1, 2 or 3")
		
	ret = .C ("l1median_CG", PACKAGE="pcaPP", NAOK = TRUE, 
		par = as.integer (c(dim (X), maxit, trace, type)),
		npar.out = integer (4),
		dpar = as.double (c(-Inf, tol)),
		dpar.out = double (1),
		as.double (X),
		#as.double (pscale),
		med = as.double (m.init)#double (ncol(X))
		)

	if (trace >= 1)
		cat ("l1median returned", ret$npa, ret$dpar, "\r\n") ;

	return (list (par = ret$med, value = ret$dpar.out[1], code = ret$npar.out [1], iterations = ret$npar.out [2], iterations_gr = ret$npar.out [3], time = ret$npar.out[4]))
}


l1median_HoCr <-
function (X, maxit = 200, tol = 10^-8, zero.tol = 1e-15, trace = 0, m.init = .colMedians (X), ...)
{
	X <- .Conv2Matrix (X, substitute (X))

	if (length (m.init) != ncol (X))
		stop (paste ("length of vector m.init (=", length (m.init), ") does not match the number of columns of data object X (=", ncol (X),")", sep = ""))

	ret.C = .C ("l1median_HoCr", PACKAGE="pcaPP"
		, npar = as.integer (c(dim (X), maxit, trace))
		, npar.out = integer (4)
		, as.double (c (tol, zero.tol))
		, as.double (X)
		, med = as.double (m.init)
		)

	if (trace >= 1)
	{
		if (ret.C$npar.out[1] == 1)
			cat ("Algorithm did not converge (return code 1).\n")
		else if (ret.C$npar.out[1] == 2)
			cat ("Step halving failed (return code 2).\n")
		else if (ret.C$npar.out[1] == 3)
			cat ("A concentration of more than n/2 observations in one point has been detected (return code 3).\n")
	}


	return (list (par = ret.C$med, value = 	sum (sqrt (colSums ((t(X) - ret.C$med)^2))),
		 code = ret.C$npar.out[1], iterations = ret.C$npar.out [2] + 1, time = ret.C$npar.out[3]))
}

l1median_NLM <-
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = .colMedians (X), ...)
{
	msg = 8
	X <- .Conv2Matrix (X, substitute (X))

	if (length (m.init) != ncol (X))
		stop (paste ("length of vector m.init (=", length (m.init), ") does not match the number of columns of data object X (=", ncol (X),")", sep = ""))

	ret = .C ("l1median_NLM", PACKAGE="pcaPP", NAOK = TRUE, 
		npar = as.integer (c(dim (X), maxit, 0, 0, 0, msg, trace)),
		dpar = as.double (c(tol, 0)),
		as.double (X),
		med = as.double (m.init) 
		)

	if (trace >= 1)
		cat ("l1median returned", ret$npa, ret$dpar, "\r\n") ;
		
		if (ret$npar[7])
			stop (paste ("nlm optimization returned error code", ret$npar[7]))

	return (list (par = ret$med, value = ret$dpar[2], code = ret$npar[4], iterations = ret$npar [3], time = ret$npar[6]))
}


.l1median_NLM_Hess <-
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = .colMedians (X), msg = 8, method = 1, GFlag = 1, HFlag = 1, Exp = 1, Digits = 6, ...)
{
	X <- .Conv2Matrix (X, substitute (X))

	if (length (m.init) != ncol (X))
		stop (paste ("length of vector m.init (=", length (m.init), ") does not match the number of columns of data object X (=", ncol (X),")", sep = ""))

	ret = .C ("l1median_NLM_Hess", PACKAGE="pcaPP", NAOK = TRUE, 
		npar = as.integer (c(dim (X), maxit, 0, method, 0, msg, trace, GFlag, HFlag, Exp, Digits)),
		dpar = as.double (c(tol, 0)),
		as.double (X),
		med = as.double (m.init)
		)

	if (trace >= 1)
		cat ("l1median returned", ret$npa, ret$dpar, "\r\n") ;

		if (ret$npar[7])
			stop (paste ("nlm optimization returned error code", ret$npar[7]))

	return (list (par = ret$med, value = ret$dpar[2], code = ret$npar[4], iterations = ret$npar [3], time = ret$npar[6]))
}


l1median_NM <-
#function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = .colMedians (X), alpha = 1, beta = 0.5, gamma = 2, ...)
function (X, maxit = 200, tol = 10^-8, trace = 0, m.init = .colMedians (X), ...)
{
	alpha = 1
	beta = 0.5
	gamma = 2
	
	X <- .Conv2Matrix (X, substitute (X))
	
	if (length (m.init) != ncol (X))
		stop (paste ("length of vector m.init (=", length (m.init), ") does not match the number of columns of data object X (=", ncol (X),")", sep = ""))
		
	ret = .C ("l1median_NM", PACKAGE="pcaPP", NAOK = TRUE, 
		npar = as.integer (c(dim (X), maxit, 0, 0, 0, 0, trace)),
		dpar = as.double (c(-Inf, tol, 0, alpha, beta, gamma)),
		as.double (X),
		#as.double (pscale),
		med = as.double (m.init)#double (ncol(X))
		)

	if (trace >= 1)
		cat ("l1median returned", ret$npa, ret$dpar, "\r\n") ;

	return (list (par = ret$med, value = ret$dpar[3], code = ret$npar[4], iterations = ret$npar [6], time = ret$npar[7]))
}

l1median_VaZh <- 
function (X, maxit = 200, tol = 10^-8, zero.tol = 1e-15, trace = 0, m.init = .colMedians (X), ...)
{

	X <- .Conv2Matrix (X, substitute (X))

	if (length (m.init) != ncol (X))
		stop (paste ("length of vector m.init (=", length (m.init), ") does not match the number of columns of data object X (=", ncol (X),")", sep = ""))

	ret.C = .C ("l1Median_VZ", PACKAGE="pcaPP"
		, npar = as.integer (c(dim (X), maxit, 0, trace))
		, nParOut = integer (3)
		, as.double (c (tol, zero.tol))
		, as.double (X)
		, med = as.double (m.init)
		)

	return (list (par = ret.C$med, value = sum (sqrt (colSums ((t(X) - ret.C$med)^2))),
		 code = ret.C$nParOut[1], iterations = ret.C$nParOut [2], time = ret.C$nParOut[3]))
}
