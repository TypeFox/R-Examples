
opt.TPO <- function (x, k.max = ncol (x), n.lambda = 30, lambda.max, ...)
{
	store.opt = TRUE
	ret <- .sPCAgrid.opt.ind (x = x, k.max = k.max, n.lambda = n.lambda, lambda.max = lambda.max, store.PCs = store.opt, f.eval = .TPO, ...)
	class (ret) <- c (class (ret), "opt.TPO")
	return (ret)
}

opt.BIC <- function (x, k.max = ncol (x), n.lambda = 30, lambda.max, ...)
{
	store.opt = TRUE
	ret <- .sPCAgrid.opt.tot (x = x, k.max = k.max, n.lambda = n.lambda, lambda.max = lambda.max, store.PCs = store.opt, f.eval = .BIC.RSS, ...)
	class (ret) <- c (class (ret), "opt.BIC")
	return (ret)
}

.flexapply <- function (X, f, NAME, args)
{
	args[[NAME]] <- X
	do.call (f, args)
}

.sPCAgrid.ml <- function (..., lambda, f.pca = .sPCAgrid.ini, f.apply = lapply)
{
	args <- list (...)
	args$store.call <- FALSE
	f.apply (X = lambda, FUN = .flexapply, f = f.pca, NAME = "lambda", args = args)
}

#.listset <- function (x, y, name = "y")
#{
#	stopifnot (length (x) == length (y))
#
#	for (i in 1:length (x))
#		x[[i]][[name]] <- y[[i]]
#	x
# }

.sPCAgrid.opt.tot <- function (x, n.lambda = 101, k.max = 2, lambda, lambda.ini, lambda.max, trace = 0, store.PCs = TRUE, f.apply = lapply, f.eval = .TPO, ...)
{
	pc.ini	<- NULL
	f.pca <- .sPCAgrid.ini
#	f.eval <- .TPO

	if (!missing (lambda.ini) && !is.null (lambda.ini))			##	use lambda.ini
	{
		k.ini <- length (lambda.ini)
		pc.ini <- f.pca (x = x, lambda = lambda.ini, k = k.ini, cut.pc = FALSE, ...)

		if (!missing (pc.ini) && !is.null (pc.ini))
			warning ("argumens \x22pc.ini\x22 AND \x22lambda.ini\x22 were specified. Ignoring \x22pc.ini\x22.")
	}
	else if (!missing (pc.ini) && !is.null (pc.ini))			##	use pc.ini
	{
		k.ini = pc.ini$k
		lambda.ini <- rep (NA, k.ini)
	}
	else 
	{
		pc.ini <- NULL
		lambda.ini <- NULL
		k.ini <- 0
	}

#	if (missing (lambda) && !missing (lambda.max))								##	Change 20120407
#		lambda.max <- rep (lambda.max, len = k.max)								##	Change 20120407

	p <- ncol (x)
	if (k.ini == p)
		stop ("all components have already been computed")

	if (k.ini + k.max > p)
	{
		warning (paste ("reducing k.max to", p - k.ini))
		k.max <- p - k.ini
	}

	if (missing (lambda))
	{
#		if (missing (lambda.max) || is.na (lambda.max[i]))						##	Change 20120407
		if (missing (lambda.max) || is.na (lambda.max))
			max.fs <- .FSgetLambda (x = x, k = k.max, pc.ini = pc.ini, f.pca = f.pca, scores = FALSE, trace = trace, ...)
		else
#			max.fs <- lambda.max[i]												##	Change 20120407
			max.fs <- lambda.max[1]

		lambda <- seq (0, max.fs, len = n.lambda)
	}

	PCs <- .sPCAgrid.ml (x = x, pc.ini = pc.ini, k = k.max, scores = FALSE, cut.pc = TRUE, trace = trace, f.pca = f.pca, lambda = lambda, f.apply = f.apply, ...)

	opt <- .SPCAgrid.opt (x, PCs, f.eval, store.PCs, k = k.ini + 1:k.max, singlePC = FALSE, ...)

	if (!store.PCs)
		return (opt)

	ret <- list (pc = list (), pc.noord = list (), 
					x = x, k.ini = k.ini, opt = opt)

	for (i in 1:length (opt$pc))
	{
		pc <- .cut.pc (opt$pc[[i]], k.ini + k.max)
		ret$pc.noord[[i]] <- pc
		ret$pc[[i]] <- .orderPCs (pc, k.max, k.ini, TRUE)
	}

	class (ret) <- "sPCAgrid.opt.tot"
	return (ret)
}

.sPCAgrid.opt.eval <- function (x, f.eval = .BIC.RSS, k = 1, ...)
{
	if (all (class (x) != "sPCAgrid.opt.tot"))
		stop ("x must be of type \"sPCAgrid.opt.tot\"")
	##2do:  check parameter k

	.SPCAgrid.opt (x$x, x$opt$PCs, f.eval, storePCs = FALSE, k = 1:k, ...)
}

.SPCAgrid.opt <- function (x, PCs, k, f.eval, store.PCs, f.apply = sapply, singlePC = TRUE, ...)
{
	ret <- list ()

	if (store.PCs)
		ret$PCs <- PCs

	if (missing (f.eval) || !is.function (f.eval))
	{
		if (!store.PCs)
			stop ("either store.PCs must be TRUE, or f.eval must be a valid model evaluation function")
		return (ret)
	}

	ret$pc <- ret$k <- list ()

	ret$mode <- .GetFunctionName (f.eval, ...)

	for (i in 1:length (k))
	{
		if (singlePC)
			K <- k[i]
		else
			K <- k[1:i]

		obj.pc.0 <- f.eval (x = x, pc = PCs[[1]], k = K, ...)
		#.l0sparse (x = x, px = PCs[[1]], k = K, ...)
		obj.pc.1 <- f.eval (x = x, pc = PCs[[length (PCs)]], k = K, obj.pc.0 = obj.pc.0, ...)
		#.sumVar (x = x, pc = PCs[[length (PCs)]], k = K, ...)

		obj <- f.apply (X = PCs, FUN = .flexapply, f = f.eval, NAME = "pc", args = list (x = x, k = K, obj.pc.0 = obj.pc.0, obj.pc.1 = obj.pc.1, ...))
		ret$obj <- cbind (ret$obj, obj)

		idx.best <- which.min (obj)
		ret$idx.best <- cbind (ret$idx.best, idx.best)
		ret$pc[[i]] <- PCs[[idx.best]]
		ret$k[[i]] <- K
	}

	if (!store.PCs)
		return (ret$pc)

	return (ret)
}

.sPCAgrid.opt.ind <- function (x, n.lambda = 101, k.max = ncol (x), lambda.ini, lambda.max, trace = 0, store.PCs = TRUE, f.eval = .TPO, ...)
{
	pc.ini	<- NULL
	f.pca <- .sPCAgrid.ini

	if (!missing (lambda.ini) && !is.null (lambda.ini))		##	go for lambda.ini
	{
		k.ini <- length (lambda.ini)
		pc.ini <- f.pca (x = x, lambda = lambda.ini, k = k.ini, cut.pc = FALSE, ...)
		if (!missing (pc.ini))
			stop ("either store.PCs must be TRUE, or f.eval must be a valid model evaluation function")
	}
	else if (!missing (pc.ini) && !is.null (pc.ini))		##	go for pc.ini
	{
		k.ini = pc.ini$k
		lambda.ini <- rep (NA, k.ini)
	}
	else 
	{
		pc.ini <- NULL
		lambda.ini <- NULL
		k.ini <- 0
	}

	if (!missing (lambda.max))
		lambda.max <- rep (lambda.max, len = k.max)

	p <- ncol (x)
	if (k.ini == p)
		stop ("all components have already been computed")

	if (k.ini + k.max > p)
	{
		warning (paste ("reducing k.max to", p - k.ini))
		k.max <- p - k.ini
	}

	if (store.PCs)
		opt <- list ()

	for (i in 1:k.max)
	{
		if (missing (lambda.max) || is.na (lambda.max[i]))
			max.fs <- .FSgetLambda (x = x, k = 1, pc.ini = pc.ini, f.pca = f.pca, scores = FALSE, trace = trace, ...)
		else
			max.fs <- lambda.max[i]

		lambda <- seq (0, max.fs, len = n.lambda)

		cur.pcs <- .sPCAgrid.ml (x = x, pc.ini = pc.ini, k = 1, scores = FALSE, cut.pc = FALSE, trace = trace, lambda = lambda, f.pca = f.pca, ...)

		cur.opt <- .SPCAgrid.opt (x, cur.pcs, f.eval, store.PCs, k = k.ini + i, ...)

		if (store.PCs)
			opt[[i]] <- cur.opt

		pc.ini <- cur.opt$pc[[1]]
	}

	pc.ini <- .cut.pc (pc.ini, k.ini + k.max)
	pc <- .orderPCs (pc.ini, k.max, k.ini, TRUE)

	if (!store.PCs)
		return (pc)

	ret <- list (pc = pc, pc.noord = pc.ini, x = x, k.ini = k.ini, opt = opt)

	class (ret) <- "sPCAgrid.opt.ind"

	return (ret)
}

.GFSL.calc.sPCA <- function (k.check, lambda = 1, f.pca = sPCAgrid, zero.tol = 1e-10, trace = 0, check.all = TRUE, ...)
{
	if (trace >= 5)
		.flush.cat ("checking lambda: ", lambda, "\r\n", sep = "")

	if (check.all)
		return (f.pca (lambda = lambda, k = k.check, ...)$loadings[, 1:k.check])

	return (f.pca (lambda = lambda, ...)$loadings[, k.check, drop = FALSE])
}

.GFSL.is.sparse <- function (k.check, zero.tol = 1e-10, ...)
{
	load <- .GFSL.calc.sPCA (k.check = k.check, zero.tol = zero.tol, ...)
	return ((sum (abs (load)> zero.tol) ) == ncol (load))
}

.GFSL.is.same <- function (sparse.load, zero.tol = 1e-10, ...)
{
	load <- .GFSL.calc.sPCA (zero.tol = zero.tol, ...)

	return (sum (abs (sparse.load - load)) <= zero.tol)
}

.GFSL.find.max <- function (lambda = 1, ...)
{
	for (i in 1:16)
	{
		if (.GFSL.is.sparse (lambda = lambda, ...))
			return (lambda)

		lambda = lambda * 2
	}
	return (NULL)	
}

.GFSL.find.range <- function (lbL , ubL, niter = 6, f.sparse = .GFSL.is.sparse, ...)
{
	for (i in 1:niter)
	{
		mbL <- mean (c(ubL, lbL))
		if (f.sparse (lambda = mbL, ...))
			ubL <- mbL
		else
			lbL <- mbL
	}
	return (ubL)
}

.getFullSparseLambda.all <- function (uBound, ...)
{	
	if (missing (uBound))
		uBound <- .GFSL.find.max (lambda = 1, check.all = TRUE, ...)

	if (is.null (uBound))
		stop ("cannot find full sparse model")

	lBound <- ifelse (uBound > 1, uBound / 2, 0)

	.GFSL.find.range (lbL = lBound, ubL = uBound, f.sparse = .GFSL.is.sparse, ...)
}

.getFullSparseLambda.indiv <- function (niter = 15, ...)
{
	uBound <- .GFSL.find.max (lambda = 1, check.all = FALSE, ...)

	if (!is.null (uBound))
	{
#		return (.getFullSparseLambda.all (uBound = uBound, ...))
		lBound <- ifelse (uBound > 1, uBound / 2, 0)

		return (.GFSL.find.range (lbL = lBound, ubL = uBound, f.sparse = .GFSL.is.sparse, check.all = FALSE, ...))
	}

	lambda.max <- 1e6

	load <- .GFSL.calc.sPCA (lambda = lambda.max, check.all = FALSE, ...)

	.GFSL.find.range (lbL = 0, ubL = lambda.max, niter = 25, f.sparse = .GFSL.is.same, sparse.load = load, check.all = FALSE, ...)
}



.FSgetLambda <- function (...)
{							##	find lambda which yields full sparseness

							##	Problem: if due to pc.ini the full sparse loadings matrix (vector)
							##		is not a 0-1 matrix (vector):
							##		In such a case full sparseness is achvieved, if a change of lambda 
							##		wouldn't change the loadings matrix (vector) any more
	kc <- .FSgetK.check (...)
	if (.FSpossible (...))
		.FSgetLambdaFS (..., k.check = kc)
	else
		.FSgetLambdaC (..., k.check = kc)
}

.FSgetK.check <- function (pc.ini = NULL, k.ini, k, ...)
{
	if (is.null (pc.ini))
		return (1:k)

	if (missing (k.ini))
		k.ini <- pc.ini$k

	return (k.ini + 1 : k)
}

.FSpossible <- function (pc.ini = NULL, k.ini, k, zero.tol = 1e-16, ...)
{
	if (is.null (pc.ini))
		return (TRUE)

	if (missing (k.ini))
		k.ini <- pc.ini$k

	l0 <- abs (pc.ini$loadings[, 1:k.ini, drop = FALSE]) > zero.tol
							##	return whether at least one row of lambda does only have zeros.
	return (any (rowSums (l0) == 0))
							
}

.FSgetLambdaFS <-function (iter = 15, ...) 
{
	.FSiter (..., f.test = .FSisFullSparse)
}

.FSgetLambdaC <- function (zero.tol, ...)
{
	pc0 <- .FScalc (..., lambda = 2^20)[[1]]
	.FSiter (..., pc0 = pc0, f.test = .FScompL)
}

.FSiter <- function (niter = 8, testL = testE^-2, testU = testE^5, testE = 16, ...)
{
	L <- testL
	U <- testU

	for (i in 1:2)
	{
		lambda.test <- testE^seq (log (L) / log (testE), log (U) / log (testE), len = niter)
		sparse <- .FStest (..., lambda = lambda.test)

		if (sparse[1])
			return (lambda.test[1])
		if (!sparse[niter])
			return (lambda.test[niter])
		idx <- which (sparse)[1]

		L <- lambda.test[idx-1]
		U <- lambda.test[idx]
	}
	return (lambda.test[idx])
}

.FStest <- function (..., f.test)
{
	##	2do. check, wether f.apply is specified. if not, do a binary search!
	PCs <- .FScalc (...)
	sapply (PCs, f.test, ...)
}

.FScompL <- function (pc, pc0, k.check, zero.tol = 1e-16, ...)
{
	sum (abs (.loadSgnU (pc$load [, k.check, drop = FALSE]) - .loadSgnU (pc0$load [, k.check, drop = FALSE])) > sqrt (zero.tol)) == 0
}

.FScalc <- function (...)	{	.sPCAgrid.ml (...)	}

.FSisFullSparse <- function (pc, k.check, zero.tol = 1e-16, ...)
{
	n.k <- length(k.check)									##	number of loadings vectors to check
	sum (abs (pc$load [, k.check]) > sqrt (zero.tol)) == n.k		##	the number of non-zero loadings is equal to n.k
}

.loadSgnU <- function (x)
{														##	change the signs of the columns of a loadings matrix such, that each column's absolute maximum is positive
	idx.max <- apply (abs (x), 2, which.max)
	sgn <- sign (x[cbind (idx.max, 1:ncol (x))])
	if (length (sgn) == 1)
		return (x * sgn)
	return (x %*% diag (sgn))
}
