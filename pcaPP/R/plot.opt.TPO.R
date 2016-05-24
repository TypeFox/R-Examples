
plot.opt.TPO <- function (x, k, f.x = c ("l0", "pl0", "l1", "pl1", "lambda"), f.y = c ("var", "pvar"), ...)
{
	.plot.sPCAgrid.opt.ind (x = x, k = k, f.x = f.x, f.y = f.y, ...)
}

plot.opt.BIC <- function (x, k, f.x = c ("l0", "pl0", "l1", "pl1", "lambda"), f.y = c ("var", "pvar"), ...)
{
	.plot.sPCAgrid.opt.tot (x = x, k = k, f.x = f.x, f.y = f.y, ...)
}

.eval.fx <- function (f.x)
{
	if (is.function (f.x))
		return (f.x)
	f.x <- f.x[[1]]

	.valf <- c ("l0", "l1", "lambda", "pl0", "pl1")

	f.x <- match.arg (f.x, .valf)
	f.x <- match (f.x, .valf)

	.F <- c (.l0sparse, .l1sparse, .lambda, .pl0sparse, .pl1sparse)

	return (.F[[f.x]])
}

.eval.fy <- function (f.y) #, supr)
{
	if (is.function (f.y))
		return (f.y)
	f.y <- f.y[[1]]

	.valf <- c("var", "dvar", "pvar", "dpvar") #, "obj") 

	f.y <- match.arg (f.y, .valf)
#	if (!missing (supr))
#		if (any (supr == f.y))
#			stop (paste ("Option \"", f.y, "\"cannot be used in this context", sep = ""))
	f.y <- match (f.y, .valf)

	.F <- c (.sumVar, .sumVarDiff, .sumVarP, .sumVarDiffP) #, .obj)

	return (.F[[f.y]])
}


.format.PCidx <- function (k, kto)
{
	if (length (k) > 1)
	{
		if (!missing (kto))
			stop ("if kto is specified, k must be of length 1")
		kto <- max (k)
		k <- min (k)
	}
	else
	{
		if (missing (kto))
			kto <- k
		else if (kto < k)
		{
			temp <- k
			k <- kto
			kto <- temp
		}
	}

	if (k == kto)
		return (paste ("PC", k))
	return (paste ("PCs ", k, "-", kto, sep = ""))
}

.title.trdoffC <- function (k, kto)
{
	paste ("Tradeoff Curve (", .format.PCidx (k, kto), ")", sep = "")
}

.plot.sPCAgrid.opt.tot <- function (x, k, main, ...) #, f.x = "l0", f.y = "pvar", ylim, xlab, ylab, main.pre, main.suf)
{
	if (missing (k))
		k = x$k.ini + length (x$pc)

	opt.idx <- k - x$k.ini

	stopifnot (length (k) == 1)
	stopifnot (any (class (x) == "sPCAgrid.opt.tot"))

	if (is.null(x$opt))
		stop ("x does not contain iteration data. Set \x22store.PCs = TRUE\x22 for function sPCAgrid.opt.tot")

	if (missing (main))
		main <- .title.trdoffC (x$k.ini + 1, k)

	ret <- .plot.opt (x = x$x, opt = x$opt, opt.idx = opt.idx, main = main, ...)
#	.opt.summary (x = x$x, pc = x$pc.noord[[opt.idx]], opt = x$opt, k = x$opt$k[[opt.idx]], ...)
	.opt.summary (x = x$x, opt = x$opt, opt.idx = opt.idx, idx.lambda = FALSE, ...)

	return (invisible (ret))
}

.plot.sPCAgrid.opt.ind <- function (x, k, main, ...) #, f.x = "l0", f.y = "pvar", ylim, xlab, ylab, main, main.pre, main.suf, ...)
{
	if (missing (k))
		k <- x$k.ini + 1

	stopifnot (length (k) == 1)
	stopifnot (any (class (x) == "sPCAgrid.opt.ind"))
	if (is.null(x$opt))
		stop ("x does not contain iteration data. Set \x22store.PCs = TRUE\x22 for function sPCAgrid.opt.ind")

	n.k <- length (x$pc$lambda)
	.assureRange (k, x$k.ini + 1, x$k.ini + n.k)

	opt.idx <- k - x$k.ini

	if (missing (main))
		main <- .title.trdoffC (k)

	opt <- x$opt[[opt.idx]]

	ret <- .plot.opt (x = x$x, opt = opt, main = main, ...)

	.opt.summary (x = x$x, pc = x$pc.noord, opt = opt, k = k, ...)

	return (invisible (ret))
}

.applyPFunc <- function (PCs, pc, idx.pc, f, ...)
{
	if (missing (pc))
		pc <- PCs[[idx.pc]]

	f.0 <- f (pc = PCs[[1]], ...)
	f.1 <- f (pc = PCs[[length (PCs)]], v.0 = f.0, ...)

	f (pc = pc, v.0 = f.0, v.1 = f.1, ...)
}

.plot.opt <- function (x, opt, opt.idx = 1, f.x = "pl0", f.y = "pvar", ylim, xlab, ylab, main, main.pre, main.suf, ...)
{
	k <- opt$k[[opt.idx]]

	f.x <- .eval.fx (f.x)
	f.y <- .eval.fy (f.y)

	if (missing (xlab))
		xlab <- .GetFunctionName (f.x, k = k)
	if (missing (ylab))
		ylab <- .GetFunctionName (f.y, k = k)

	PCs <- opt$PCs

	pdat <- .trdoff.sPCAgrid.ind (x, PCs, k = k, f.x, f.y, ...)

	if (missing (main))
		main <- paste ("Tradeoff Curve PC", k)
	if (!missing (main.pre))
		main <- paste (main.pre, main)
	if (!missing (main.suf))
		main <- paste (main, main.suf)

	if (missing (ylim))
		ylim <- c (0, max (pdat [,2]))
	plot (pdat[order(pdat[,1]),], type = "b", xlab= xlab, ylab = ylab, ylim = ylim, main = main)

	if (is.null (opt$idx.best))
		return (invisible (pdat))

	abline (v = pdat[opt$idx.best[opt.idx], 1], lty = 2)

	return (invisible (pdat))
}

.opt.summary <- function (x, pc, optmode, k, opt, opt.idx, idx.lambda = TRUE, ...)
{
	if (missing (pc))
		pc <- opt$pc[[opt.idx]]
	if (missing (k))
		k <- opt$k[[opt.idx]]
	if (missing (optmode))
		optmode <- opt$mode

	lam <- round (pc$lambda [max (k)], 2)

	if (length (k) > 1)
		k.use <- k
	else
		k.use <- 1:k

	vartot <- .sumVarP (x = x, pc = pc, k = k.use, ...)
	ecv <- .sumVarP  (x = x, pc = pc, k = k.use, v.0 = vartot, ...)
#	ecv <- .applyPFunc (opt$PCs, pc, f = .sumVarP, k = k.use, x = x, ...)
	ecv <- round (ecv, 2)
	
	LoS <- .l0sparse (pc = pc,k = k.use, ...)

#	LoS <- .applyPFunc (opt$PCs, pc, f = .l0sparse, k = k.use, x = x, ...)

	if (length (k.use) > 1)
		txtK <- paste (min(k.use), "-", max(k.use), sep = "")
	else
		txtK <- k.use

	txtECV <- bquote (paste ("ECV" [.(txtK)], ": ", .(ecv), "%"))
	txtLoS <- bquote (paste ("L"[0], "S" [.(txtK)], ": ", .(LoS)))


	if (idx.lambda)
	{
		ssL <- paste ("opt", max (k), sep = "")
		txtLambda <- bquote (paste (lambda[.(ssL)], ": " ,.(lam)))
	}
	else
		txtLambda <- bquote (paste (lambda[opt], ": " ,.(lam)))

	txt <- bquote(paste (.(optmode), " - ", .(txtLambda), "; ", .(txtECV), "; ", .(txtLoS)))

	mtext (txt, line = 0.25, cex = 0.8)
}

.trdoff.sPCAgrid.ind <- function (x, PCs, f.x, f.y, ...)
{
	n <- length (PCs)

	x.0 <- f.x (x = x, pc = PCs[[1]], ...)
	x.1 <- f.x (x = x, pc = PCs[[n]], v.0 = x.0, ...)

	y.0 <- f.y (x = x, pc = PCs[[1]], ...)
	y.1 <- f.y (x = x, pc = PCs[[n]], v.0 = y.0, ...)

	rx <- sapply (PCs, .flexapply, f = f.x, NAME = "pc", args = list (x = x, v.0 = x.0, v.1 = x.1, ...))
	ry <- sapply (PCs, .flexapply, f = f.y, NAME = "pc", args = list (x = x, v.0 = y.0, v.1 = y.1, ...))

	stopifnot (is.numeric (rx))
	stopifnot (is.numeric (ry))

	ret <- cbind (rx, ry)
#	ret <- ret [order (ret[,1]), ]

	class (ret) <- "trdoff.sPCAgrid"
	ret 
}

.assureRange <- function (x, l, u, name = substitute (x))
{
	if (missing (l))
	{
		if (x > u)
			stop (paste (name, "must not be larger than ", u))
	}
	else if (missing (u))
	{
		if (x < l)
			stop (paste (name, "must not be smaller than ", l))
	}
	else if (x < l || x > u)
		stop (paste (name, "must be between", l, "and", u))
}

objplot <- function (x, k, ...)
{
	if (any (class (x) == "sPCAgrid.opt.tot"))
	{
		if (missing (k))
			k = x$k.ini + length (x$pc)

		.assureRange (k, x$k.ini, x$k.ini + length (x$pc))

		opt.idx <- k - x$k.ini
		.objplot (x$opt, opt.idx = opt.idx)
		.opt.summary (x = x$x, opt = x$opt, opt.idx = opt.idx, idx.lambda = FALSE, ...)
	}
	else if (any (class (x) == "sPCAgrid.opt.ind"))
	{
		if (missing (k))
			k = x$k.ini + 1
		.assureRange (k, x$k.ini, x$k.ini + length (x$pc$lambda))

		opt <- x$opt[[k - x$k.ini]]
		.objplot (opt, k = k)
		.opt.summary (x = x$x, pc = x$pc.noord, opt = opt, k = k, ...)
	}
}

.objplot <- function (opt, opt.idx = 1, k, main, xlab, ylab, type = "b", ...)
{
	obj <- opt$obj[,opt.idx]

	if (missing (k))
		k <- opt$k[[opt.idx]]

	lambda <- sapply (opt$PCs, get ("[["), "lambda")
	if (is.matrix (lambda))
		lambda <- lambda [nrow (lambda),]

	if (missing (main))
		if (length (k) > 1)
			main <- paste ("Model Selection (PCs ", min (k), "-", max(k), ")", sep = "")
		else
			main <- paste ("Model Selection (PC ", k, ")", sep = "")

	if (missing (xlab))
		xlab <- expression (lambda)
	if (missing (ylab))
		ylab <- paste ("Objective Function (", opt$mode, ")", sep = "")

	plot (lambda, obj, xlab = xlab, ylab = ylab, main = main, type = type, ...)

	cl <- lambda [which.min (obj)]

	abline (v = cl, lty = 2)

#	txt <- bquote (lambda[opt] == .(round (cl, 2)))
#	mtext (line = 0.25, text = txt, cex = 0.8)

	invisible ()
}
