# Created by David Bickel on 18 July 2007.

#library(locfdr)
#try(library(limma)) # library(statmod) ?
#
#
## Source("data.r") #reprod.r
#Source("estimate.r") # for unsortableEstimate; changed from "data.r" on 6 November 2007
#Source("HighProbability1.s") # for bFdr
#Source("MinProbability.s") # for bFdr
setClass("MArrayLM", representation("Numeric"))#marta: I cannot find this definition anywhere, so I created it.
setClass("basicFdr", representation(fdr = "Numeric", p.value = "Numeric"))
setValidity("basicFdr", function(object)
{
	nam <- names(object) # as(object, "Numeric"))
	nam.ok <- !is.null(nam) && sameNames(object@fdr, object@p.value)
	range.ok <- all(c(object@fdr, object@p.value) <= 1.001, na.rm = TRUE)
	ok <- nam.ok && range.ok
	if(!ok)
	{ message("invalid basicFdr"); browser()}
	ok
})
setClassUnion("Fdr", c("basicFdr"))
setAs(from = "Fdr", to = "basicFdr", function(from)
{
	new.basicFdr(fdr = fdr(from), p.value = pvalue(from))
})
setMethod("length", signature(x = "basicFdr"), function(x)
{
	length(names(x))
})
setMethod("names", signature(x = "basicFdr"), function(x)
{
	names(x@p.value)
})
setMethod("[", signature(x = "basicFdr", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@fdr <- (x@fdr)[i]
  x@p.value <- (x@p.value)[i]
  x
})
setMethod("annotation", signature(object = "basicFdr"), function(object)
{
	"basic"
})

setClass("eFdr", representation("list", fdr = "Numeric", zz = "numeric")) # e is for Efron
setValidity("eFdr", function(object)
{
	nam <- names(object) # as(object, "Numeric"))
	nam.ok <- !is.null(nam) && sameNames(object@fdr, object@zz)
#	if(nam.ok)
#		message("Some names: ", paste(nam[1:10], collapse = "; "))
	ok <- nam.ok
	if(!ok)
	{ message("invalid eFdr"); browser()}
	ok
})
setMethod("length", signature(x = "eFdr"), function(x)
{
	length(names(x))
})
setMethod("names", signature(x = "eFdr"), function(x)
{
	names(x@zz) # names(as(x, "Numeric"))
})
setMethod("[", signature(x = "eFdr", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@fdr <- (x@fdr)[i]
  x@zz <- (x@zz)[i] # error fixed 27 October 2007
  x
})
setMethod("plot", signature(x = "eFdr", y = "missing"), function(x, y, call.par, ...) # added 12 Feb. 2008
{
	z <- x@zz
	plot.zz <- function(fdr.lab, log)
	{
		plot(x = z, y = fdr(x), xlab = "z", ylab = fdr.lab, log = log)
	}
	function.list <- list(plot.zz)
	plot(x = as(x, "basicFdr"), function.list = function.list, ...)
})
setMethod("annotation", signature(object = "eFdr"), function(object)
{
	"locfdr"
})

# bFdr added 15 Feb. 2008:
setClass("bFdr", representation(fdr = "Numeric", p.value = "Numeric", prior.fdr = "Scalar")) # b is for Bickel
setValidity("bFdr", function(object)
{
	nam <- names(object) # as(object, "Numeric"))
	nam.ok <- !is.null(nam) && sameNames(object@fdr, object@p.value)
	range.ok <- all(c(object@fdr, object@p.value, object@prior.fdr) <= 1.001, na.rm = TRUE)
	ok <- nam.ok && range.ok
	if(!ok)
	{ message("invalid bFdr"); browser()}
	ok
})
setMethod("length", signature(x = "bFdr"), function(x)
{
	length(names(x))
})
setMethod("names", signature(x = "bFdr"), function(x)
{
	names(x@p.value) # names(as(x, "Numeric"))
})
setMethod("[", signature(x = "bFdr", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@fdr <- (x@fdr)[i]
  x@p.value <- (x@p.value)[i]
  x
})
setMethod("annotation", signature(object = "bFdr"), function(object)
{
	"empiricalBayes"
})

setClass("pseudoFdr", representation(fdr = "Numeric", estimate = "numeric", FUN = "function"))
setValidity("pseudoFdr", function(object)
{
	nam <- names(object)
	nam.ok <- !is.null(nam) && sameNames(object@fdr, object@estimate)
	ok <- nam.ok
	if(!ok)
	{ message("invalid pseudoFdr"); browser()}
	ok
})
setMethod("length", signature(x = "pseudoFdr"), function(x)
{
	length(names(x))
})
setMethod("names", signature(x = "pseudoFdr"), function(x)
{
	names(x@estimate)
})
setMethod("[", signature(x = "pseudoFdr", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@fdr <- (x@fdr)[i]
  x@estimate <- (x@estimate)[i]
  x
})
setMethod("annotation", signature(object = "pseudoFdr"), function(object)
{
	"fold-change"
})

setClass("sFdr", representation("list", fdr = "Numeric", p.value = "Numeric", t = "numeric", proportion = "Scalar")) # s is for Smyth
setValidity("sFdr", function(object)
{
	nam <- names(object) # as(object, "Numeric"))
	nam.ok <- !is.null(nam) && sameNames(object@fdr, object@p.value, object@t)
#	if(nam.ok)
#		message("Some names: ", paste(nam[1:10], collapse = "; "))
	ok <- nam.ok
	if(!ok)
	{ message("invalid sFdr"); browser()}
	ok
})
setMethod("length", signature(x = "sFdr"), function(x)
{
	length(names(x))
})
setMethod("names", signature(x = "sFdr"), function(x)
{
	names(x@p.value) # names(as(x, "Numeric"))
})
setMethod("[", signature(x = "sFdr", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@fdr <- (x@fdr)[i]
  x@p.value <- (x@p.value)[i]
  x@t <- (x@t)[i]
  x
})
setMethod("annotation", signature(object = "sFdr"), function(object)
{
	"limma"
})

setClassUnion("Fdr", c("sFdr", "eFdr", "bFdr", "pseudoFdr", "basicFdr"))
setMethod("plot", signature(x = "Fdr", y = "missing"), function(x, y, call.par, log, function.list, pvalue.lab, ...)
{
	if(missing(pvalue.lab))
		pvalue.lab <- "p-value"
	if(missing(function.list))
		function.list <- list()
	stopifnot(is.list(function.list))
	fdr.lab <- "local false discovery rate"
	if(missing(call.par))
		call.par <- TRUE
	if(call.par)
		par(mfrow = c(2, 2))
	if(missing(log))
	{
		subplot <- function(log){plot(x = x, call.par = FALSE, log = log, function.list = function.list, pvalue.lab = pvalue.lab, ...)}
		subplot(log = "")
		subplot(log = "y")
	}
	else
	{
		plot(x = pvalue(x), y = fdr(x), xlab = pvalue.lab, ylab = fdr.lab, log = log, ...)
		if(length(function.list) >= 1)
		{
			for(i in length(function.list))
				function.list[[i]](fdr.lab = fdr.lab, log = log)
		}
	}
})
setMethod("plot", signature(x = "Fdr", y = "Fdr"), function(x, y, xlab, ylab, ...)
{
	get.lab <- function(object, lab)
	{
		if(is.null(lab))
			lab <- annotation(object)
		paste("local false discovery rate (", lab, ")", sep = "")
	}
	xlab <- get.lab(object = x, lab = if(missing(xlab)) NULL else xlab)
	ylab <- get.lab(object = y, lab = if(missing(ylab)) NULL else ylab)
	plot(x = fdr(x), y = fdr(y), xlab = xlab, ylab = ylab, ...)
	abline(a = 0, b = 1, col = "gray")
})

setClass("Prob0", representation("Numeric", Scale0 = "Numeric", Scale1 = "Numeric", Location1 = "numeric", xSize = "Numeric", ySize = "Numeric", PValue = "Numeric", mle = "logical", annotation = "character"))
setValidity("Prob0", function(object)
{
	sn <- setdiff(slotNames(object), c("mle", ".Data", "annotation"))
	if(length(object@PValue) == 0)
		sn <- setdiff(sn, c("PValue"))
	len.ok <- all(length(object) == sapply(sn, function(name){length(slot(object = object, name = name))})) && length(object@mle) == 1 && length(object@annotation) <= 1
	nam.ok <- len.ok && all(sapply(sn, function(name){all(names(object) == names(slot(object = object, name = name)))}))
	val.ok <- all(object <= 1.001, na.rm = TRUE)
	ok <- nam.ok && val.ok
	if(!ok){printInvalid(object); browser()}
	ok
})
setMethod("[", signature(x = "Prob0", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@.Data <- as(x, "Numeric")[i]
  x@Scale0 <- x@Scale0[i]
  x@Scale1 <- x@Scale1[i]
  x@Location1 <- x@Location1[i]
  x@xSize <- x@xSize[i]
  x@ySize <- x@ySize[i]
  x@PValue <- x@PValue[i]
#	stop('"[", signature(x = "Prob0", i = "ANY", j = "missing") not yet implemented')
  x
})
setMethod("names", signature(x = "Prob0"), function(x)
{
	names(x@Location1)
})
setAs(from = "Prob0", to = "numeric", function(from)
{
	vec <- from@.Data
	names(vec) <- names(from)
	vec
})
setAs(from = "Prob0", to = "Numeric", function(from)
{
	vec <- as(from, "numeric")
	Numeric(vec)
})
setAs(from = "Prob0", to = "Fdr", function(from) # new 24 April 2008
{
	Num <- from@.Data
	if(is(Num, "Fdr"))
		Num
	else
	{
		warning("Prob0 may have lost Fdr info. such as eFdr's zz slot")
		new.basicFdr(fdr = Num, p.value = from@PValue)
	}
})
setMethod("plot", signature(x = "missing", y = "Prob0"), function(x, y, call.par, plot.Scale, legend.x, col, qRatio, lower.tail, log, ...)
{
	x <- stop('"plot", signature(x = "missing", y = "Prob0") not yet implemented')
#	subplot <- stop("")
})
setMethod("plot", signature(x = "Prob0", y = "biasEstimate"), function(x, y, ...)
{
	plot(x = as(x, "numeric"), y = y, xlab = "empirical probability of equivalent expression", ...)
})
setMethod("plot", signature(x = "Prob0", y = "missing"), function(x, y, call.par, plot.Scale, legend.x, col, qRatio, lower.tail, x.lower.tail, log, ProbabilityEstimate.as.x, all.estimates, file, main, ...)
{
	if(all(x == 1, na.rm = TRUE))
	{
		message("All probabilities of the nulls equal 1, so there is nothing to plot")
		return()
	}
	if(missing(main))
	{
		main <- paste(annotation(x), Sys.Date(), sep = "; ")
		message("main = ", main)
	}
	if(!missing(file))
	{
		message("argument file was supplied, so arguments other than qRatio, log, and all.estimates will be ignored")
		if(missing(log)) log <- character(0)
		if(missing(all.estimates)) all.estimates <- FALSE
		if(missing(qRatio)) stop("qRatio must be supplied when file is supplied in argument list")
		if(length(qRatio) != 4) warning("for better results with file argument specified, use qRatio of length 4")
		pdf(file = file)
		message("plotting qRatio:")
		print(plot(x = x, log = log, qRatio = qRatio, all.estimates = all.estimates, lower.tail = c(TRUE, FALSE)), main = main)
		par(mfrow = c(2, 2))
		for(qr in qRatio)
		{
			message("\nplotting qr:", qr)
			plot(x = x, log = log, ProbabilityEstimate.as.x = TRUE, qRatio = qr, all.estimates = all.estimates, main = main)
		}
		dev.off()
		return(qRatio)
	}
	advise.other.method <- if(missing(ProbabilityEstimate.as.x))
	{
		ProbabilityEstimate.as.x <- FALSE
		message("ProbabilityEstimate.as.x = ", ProbabilityEstimate.as.x)
		TRUE
	}
	else
		ProbabilityEstimate.as.x
	if(missing(x.lower.tail))
	{
		x.lower.tail <- c(TRUE, FALSE) # logical(0)
		if(ProbabilityEstimate.as.x)
			message("x.lower.tail = ", if(length(x.lower.tail) == 0) "[object of length 0]" else paste(x.lower.tail, collapse = ", "))
	}
	if(missing(all.estimates))
	{
		all.estimates <- TRUE
		message('x = "Prob0", y = "missing": all.estimates = ', all.estimates)
	}
	if(missing(plot.Scale) && (missing(lower.tail) || (is.logical(lower.tail) && length(lower.tail) == 0)))
	{
		if(advise.other.method)
			message('Also consider "plot", signature(x = "missing", y = "Prob0") instead.')
		if(missing(log)) log <- character(0)
		if(missing(call.par)) call.par <- TRUE # missing(plot.Scale)
		if(call.par) par(mfrow = c(2, 2))
		if(missing(qRatio))
		{
			qRatio <- 1.4
			message("superplot qRatio = ", qRatio)
		}
		subplot <- function(log, ...)
		{
			plot(x = x, log = log, ProbabilityEstimate.as.x = ProbabilityEstimate.as.x, x.lower.tail = x.lower.tail, all.estimates = all.estimates, main = main, ...)
		}
		tr <- subplot(log = log, qRatio = qRatio, plot.Scale = TRUE, ...)
#		if(is(tr, "try-error"))
#		{
#			message("scale plot failed")
#			y.in.log <- log == "xy" || log == "y"
#			if(y.in.log)
#			{
#				scale.log <- if(log == "xy") "x" else ""
#				message("try plotting again with log = '", log, "'")
#			}
#		}
		subplot(log = "y", qRatio = qRatio, plot.Scale = FALSE, ...)
		if(missing(qRatio))
		{
			qRatio <- 1
			message("plot qRatio = ", qRatio)
		}
		if(missing(lower.tail))
		{
			subplot(log = log, qRatio = 1 / qRatio, lower.tail = TRUE, ...)
			subplot(log = log, qRatio = qRatio, lower.tail = FALSE, ...)
		}
		else
			subplot(log = log, qRatio = qRatio, lower.tail = if(length(lower.tail) == 0) NULL else lower.tail, ...)
	}
	else # usually this part applies only if subplot is called
	{
		if(ProbabilityEstimate.as.x)
		{
			if(missing(qRatio))
			{
				qRatio <- 1.4
				message("subplot qRatio = ", qRatio)
			}
			x.qRatio <- if((qRatio >= 1 || length(x.lower.tail) == 1) && (missing(lower.tail) || length(lower.tail) != 1 || !lower.tail))
				qRatio
			else
				1 / qRatio
		#	ProbabilityEstimate.lab(qRatio = x.qRatio, lower.tail = x.lower.tail)
		}
		if(!missing(lower.tail))
		{
			stopifnot(is.logical(lower.tail) && length(lower.tail) <= 2)
			qRatio <- if(length(qRatio) == 1) Scalar(qRatio) else Numeric(qRatio) # Scalar(qRatio)
		}
		if(!missing(lower.tail) && is.null(lower.tail))
			lower.tail <- logical(0)
		if(missing(legend.x)) legend.x <- "center"
		if(missing(col))
			col <- 1:length(x)
		if(length(qRatio) > 1)
		{
			par(mfrow = c(2, 2))
			for(qr in qRatio)
			{
				message("plotting for qRatio of ", qr)
				if(missing(lower.tail))
				{ message("arg combination not yet implemented: length(qRatio) > 1 && missing(lower.tail)"); browser()}
				if(missing(log) || length(log) == 0)
					log <- ""
				plot(x = x, call.par = FALSE, plot.Scale = plot.Scale, legend.x = legend.x, col = col, lower.tail = lower.tail, x.lower.tail = x.lower.tail, ProbabilityEstimate.as.x = ProbabilityEstimate.as.x, all.estimates = all.estimates, log = log, qRatio = qr, main = main, ...)
			}
			return(qRatio)
		}
		locPE <- function(lower.tail, qRatio, ...)
		{
			lpe.ok <- qRatio >= 1 || length(lower.tail) == 1
			if(!lpe.ok)
			{ message("locPE error"); browser()}
			ProbabilityEstimate(object = x, lower.tail = lower.tail, qRatio = qRatio, ...)
		}
		get.vec <- if(missing(lower.tail))
		{
			function(...)
			{
				fun <- if(plot.Scale)
					Scale
				else
					function(...){exp(Location(...))}
				fun(object = x, ...)
			}
		}
		else
		{
			function(...)
			{
				locPE(lower.tail = lower.tail, qRatio = qRatio, ...)
			}
		}
		if(all.estimates)
		{
			vec0 <- get.vec(null.hypothesis = TRUE)
			vec1 <- get.vec(null.hypothesis = FALSE)
		}
		vec <- get.vec()
		pch0 <- 6
		pch1 <- 2
		pch <- 1
		all.y <- if(all.estimates) c(vec0, vec1, vec) else vec
		boo <- is.finite(vec)
		xvec <- if(ProbabilityEstimate.as.x)
			locPE(lower.tail = x.lower.tail, qRatio = x.qRatio, ...)
		else
			x
		stopifnot(length(xvec) == length(vec))
		xlab <- if(is(xvec, "ProbabilityEstimate"))
			annotation(xvec)
		else if(is(xvec, "Prob0"))
			"local false discovery rate"
		else
		{	message("unexpected xvec class"); browser()}
		ylab <- if(missing(lower.tail))
		{
			if(plot.Scale) "variability" else "ratio"
		}
		else if(is(vec, "ProbabilityEstimate"))
		{
			annotation(vec)
		}
		else
		{	message("unexpected vec class"); browser()}
		locplot <- function(log)
		{
			plot(x = xvec[boo], y = vec[boo], xlab = xlab, ylab = ylab, ylim = range(all.y[is.finite(all.y)], na.rm = TRUE), log = log, col = col, main = main, ...)
		}
		if(missing(log) || length(log) == 0)
			log <- if(missing(lower.tail)) "y" else ""
		tr <- try(locplot(log = log))
		if(is(tr, "try-error"))
		{ 
			message("plot ", class(x), " error")
			y.in.log <- log == "xy" || log == "y"
			if(y.in.log)
			{
				log <- if(log == "xy") "x" else ""
				message("log changed to '", log, "'")
				tr <- try(locplot(log = log))
			}
			if(!y.in.log || is(tr, "try-error"))
			{
				message("error not due to log on y-axis")
				browser()
			}
			else
				message("try plotting again with log = '", log, "'")
		}
		if(all.estimates)
		{
			add.points <- function(null.hypothesis)
			{
				y <- if(null.hypothesis) vec0 else vec1
				if(length(x) != length(y))
				{ message("cannot add points"); browser()}
				points(x = xvec, y = y, pch = if(null.hypothesis) pch0 else pch1, col = col)
			}
			add.points(TRUE)
			add.points(FALSE)
			leg <- c("reduced model estimate", "empirical Bayes estimate", "full model estimate") # c("estimate under null", "empirical Bayes estimate", "estimate under alternative")
			if(is.character(legend.x) && length(legend.x) == 1)
				legend(x = legend.x, legend = leg, pch = c(pch0, pch, pch1))
		}
	} # end subplot else
})
setMethod("annotation", signature(object = "Prob0"), function(object){object@annotation})

setClass("Probs0.biasEstimates", representation(P0 = "list", bias = "list"))
setValidity("Probs0.biasEstimates", function(object)
{
	len.ok <- length(object@P0) == length(object@bias)
	P0.nam <- names(object@P0)
	bias.nam <- names(object@bias)
	nam.ok <- len.ok && ((is.null(P0.nam) && is.null(bias.nam)) || all(P0.nam == bias.nam))
	cla.ok <- len.ok && nam.ok
	ok <- nam.ok && cla.ok
	if(!ok)
	{	printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "Probs0.biasEstimates", y = "missing"), function(x, y, file, main, ...)
{
	if(missing(main))
		main <- character(0)
	main.ok <- is.character(main) && length(main) %in% c(0, length(x@P0))
	if(!main.ok)
	{ message("main problem"); browser()}
	call.pdf <- !missing(file)
	if(call.pdf)
		stopifnot(is.character(file) && length(file) == 1)
	else
		message("specify file to plot to PDF")
	if(call.pdf)
		pdf(file)
	par(mfrow = c(2, 2))
	for(i in 1:length(x@P0))
	{
		if(i > 1 && (i - 1) %% 4 == 0 && !call.pdf)
		{
			nx11()
			par(mfrow = c(2, 2))
		}
		subplot <- function(...){plot(x = x@P0[[i]], y = x@bias[[i]], save.space = TRUE, include.estimate = FALSE, ...)}
		if(length(main) == 0)
			subplot(...)
		else
			subplot(main = main[i], ...)
	}
	if(call.pdf)
		dev.off()
})

setClass("loss", representation(lossBayes = "numeric", loss0 = "numeric", loss1 = "numeric", P0 = "Prob0", FUN = "function"))
setValidity("loss", function(object)
{
	is.name.ok <- function(x)
	{
		nam <- names(x)
		is.character(nam) && !any(is.na(nam))
	}
	loss.ok <- function(vec){is.name.ok(vec) && length(vec) == length(object@P0) && all(names(vec) == names(object@P0))}
	nam.ok <- is.name.ok(object@P0)
	l.ok <- loss.ok(object@lossBayes) && loss.ok(object@loss0) && loss.ok(object@loss1)
	nl.ok <- nam.ok && l.ok # && length(object@annotation) == 1
	if(!identical(nl.ok, TRUE))
	{ printInvalid(object); browser()}
	nl.ok
})
setAs(from = "loss", to = "predictionError", function(from)
{
	mat <- as.matrix(data.frame(from@loss0, from@lossBayes, from@loss1))
	colnames(mat) <- c("null hypotheses", "empirical Bayes", "unbiased estimation")
 	pE <- try(new("predictionError", mat, parameter = colnames(mat), parameter.lab = character(0)))
 	if(is(pE, "try-error"))
 	{ message("pE err"); browser()}
 	pE
})
setMethod("annotation", signature(object = "loss"), function(object){annotation(object@P0)})
setMethod("plot", signature(x = "loss", y = "missing"), function(x, y, main, plot.predictionError, call.par, cumulative, nfeatures, FUN, legend.x, ...)
{
	if(missing(legend.x)) legend.x <- "top"
	if(missing(nfeatures))
	{
		nfeatures <- min(length(x@P0), 500)
		message(class(x), " nfeatures = ", nfeatures)
	}
	if(missing(cumulative))
	{
		cumulative <- TRUE
		message(class(x), " cumulative = ", cumulative)
	}
	if(missing(FUN))
	{
		FUN.name <- "exp_root" # if(cumulative) "identity" else "exp_root"
		message(class(x), " FUN = ", FUN.name)
		FUN <- eval(parse(text = FUN.name))
	}
	if(missing(main)) main <- character(0)
	if(missing(plot.predictionError))
	{
		plot.predictionError <- logical(0)
		message("default plot.predictionError")
	}
	main0 <- annotation(x)
	get.main <- function()
	{
		if(length(main) == 0)
			main0
		else
			paste(main0, main, sep = "; ")
	}
	if(!is.logical(plot.predictionError) || length(plot.predictionError) == 0)
	{
		if(missing(call.par))
			call.par <- TRUE
		if(call.par)
			par(mfrow = c(2, 2))
		recurse <- function(plot.predictionError){plot(x = x, main = main, plot.predictionError = plot.predictionError, cumulative = cumulative, nfeatures = nfeatures, FUN = FUN, legend.x = legend.x, ...)}
		recurse(plot.predictionError = TRUE)
		recurse(plot.predictionError = FALSE)
	}
	else if(plot.predictionError)
	{
		message("plot sorting ", class(x))
		nsig <- sum(x@P0 < 1, na.rm = TRUE)
		if(nsig < nfeatures)
			warning(paste("Only ", nsig, " features have nonzero probabilities of differential expression, but ", nfeatures, " features are plotted.", sep = ""))
		x <- sort(x, by.PValue = FALSE)
		message("coercing ", class(x))
		x <- as(x, "predictionError")
		message("plotting ", class(x))
		plot(x = x, main = get.main(), xlab = "rank", cumulative = cumulative, nfeatures = nfeatures, FUN = FUN, use.lower.indices = TRUE, legend.x = legend.x, ...)
	}
	else
	{
		loss0 <- x@loss0
		loss1 <- x@loss1
		lossBayes <- x@lossBayes
		col1 <- "gray"
		pch1 <- 6
		colBayes <- "black"
		pchBayes <- 2
		tr <- try(plot(x = loss0, y = loss1, xlab = "null hypothesis prediction error", ylab = "prediction error", pch = pch1, col = col1, main = get.main(), log = "xy", ...))
		if(is(tr, "try-error"))
		{ message("problem plot"); browser()}
		tr.p <- try(points(x = loss0, y = lossBayes, pch = pchBayes, col = colBayes))
		if(is(tr.p, "try-error"))
		{ message("problem points"); browser()}
		abline(a = 0, b = 1, col = "blue")
		legend(x = "topleft", legend = c("unbiased estimation", "empirical Bayes"), pch = c(pch1, pchBayes), col = c(col1, colBayes))
	}
})
setMethod("plot", signature(x = "loss", y = "loss"), function(x, y, xlab, ylab, cumulative, nfeatures, indices, FUN, call.sort, main, relative, call.par, legend.x, max.level, include.individual.plots, ...)
{
	if(missing(include.individual.plots))
	{
		include.individual.plots <- TRUE
		message("include.individual.plots = ", include.individual.plots)
	}
	if(missing(max.level))
	{
		max.level <- 0.8
	}
	if(missing(legend.x)) legend.x <- "topright"
	if(missing(cumulative))
	{
		cumulative <- TRUE
		message(class(x), " ", class(y), " cumulative = ", cumulative)
	}
	if(missing(call.sort))
	{
		call.sort <- cumulative
		message("call.sort = ", call.sort)
	}
	if(missing(FUN))
	{
		FUN.name <- "exp_root" # if(cumulative) "identity" else "exp_root"
		message(class(x), " ", class(y), " FUN = ", FUN.name)
		FUN <- eval(parse(text = FUN.name))
	}
	if(missing(main)) main <- character(0)
	main0 <- if(length(annotation(x)) == 1 && annotation(x) == annotation(y))
		annotation(x)
	else
		paste(annotation(x), annotation(y), sep = " & ")
	get.main <- function(relative)
	{
		ma <- if(length(main) == 0)
			main0
		else
			paste(main0, main, sep = "; ")
		if(relative)
			ma <- paste(ma, "relative")
		ma
	}
	if(missing(indices))
	{
		if(missing(nfeatures))
		{
			nfeatures <- min(length(x@P0), length(y@P0), 500)
			message(class(x), " ", class(y), " nfeatures = ", nfeatures)
		}
		indices <- 1:nfeatures
	}
	subplot <- function(relative)
	{
		get.loss <- function(object)
		{
			get.sorted <- function(object)
			{
#				message("sorting ", class(object))
				nsig <- sum(object@P0 < 1, na.rm = TRUE)
				if(nsig < nfeatures)
					warning(paste("Only ", nsig, " features have nonzero probabilities of differential expression, but ", nfeatures, " features are plotted.", sep = ""))
				sort(object, by.PValue = FALSE, verbose = FALSE)
			}
			if(!validObject(object))
			{ message(class(object), " is invalid before any sort"); browser()}
			if(call.sort)
			{
				object0 <- object
				object <- get.sorted(object)
				message(class(object), " sorted")
				vo0 <- try(validObject(object))
				if(identical(vo0, TRUE))
					message(class(object), " sorted successfully.")
				else
				{ message("object is invalid after sort"); browser()}
			}
			acc <- accumulate(object, cumulative = cumulative, indices = indices, relative = relative)
			lo <- if(relative)
				acc@lossBayes # / (Location(acc@P0, null.hypothesis = FALSE)) ^ 2 # acc@loss0
			else
				FUN(acc@lossBayes)
			cat("loss range: "); print(range(lo, na.rm = TRUE))
			lo
		}
		xvec <- get.loss(x)
		yvec <- get.loss(y)
		boo <- is.finite(xvec) & is.finite(yvec)
		xlim <- ylim <- range(c(xvec[boo], yvec[boo]))
		col <- gray(level = max.level * (indices - 1) / length(indices))
		get.lab <- function(lab){paste(lab, if(relative) "relative error" else "error")}
		plot(x = xvec, y = yvec, main = get.main(relative = relative), xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, col = col, log = if(relative) "xy" else "", ...)
		abline(a = 0, b = 1, col = "yellow")
		stopifnot(length(col) == length(indices))
		legend(x = legend.x, legend = c(paste("rank", c(indices[1])), paste("ranks ", indices[1], "-", indices[length(indices)], sep = "")), col = c(col[1], col[length(indices)]), pch = c(1, 1))
	}
	call.two.joint.plots <- missing(relative) || !is.logical(relative) || length(relative) == 0
	if(call.two.joint.plots || include.individual.plots)
	{
		if(missing(call.par))
		{
			call.par <- TRUE
			message("relative unspecified; call.par = ", call.par)
		}
		if(call.par)
			par(mfrow = c(2, 2))
	}
	if(include.individual.plots)
	{
		individual.plot <- function(x, main){plot(x = x, main = main, plot.predictionError = TRUE)}
		individual.plot(x, main = xlab)
		individual.plot(y, main = ylab)
	}
	if(call.two.joint.plots)
	{
		subplot(relative = FALSE)
		subplot(relative = TRUE)
	}
	else
		subplot(relative = relative)
})
setMethod("[", signature(x = "loss", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@lossBayes <- x@lossBayes[i]
  x@loss0 <- x@loss0[i]
  x@loss1 <- x@loss1[i]
  x@P0 <- x@P0[i]
  x
})
setMethod("sort", signature(x = "loss"), function(x, decreasing = logical(0), verbose, by.Location1, ...)
{#XXX|:removed decreasing = "ANY"
	if(missing(verbose)) verbose <- FALSE
	if(verbose) message("\nbeginning sort ", class(x))
	if(missing(by.Location1))
	{
		by.Location1 <- FALSE
		message("sort loss by.Location1 = ", by.Location1)
	}
	if(missing(decreasing) || length(decreasing) == 0)
	{
		decreasing <- by.Location1 # FALSE
		if(verbose)
			message("decreasing = ", decreasing)
	}
	else if(verbose)
		message("sort called with actual decreasing argument of ", decreasing)
	if(decreasing != by.Location1)
	{
		message("unusual actual arguments in sort ", class(x), "; press c to continue anyway")
		browser()
	}
	ord <- Order(x, decreasing = decreasing, by.Location1 = by.Location1, ...)
	miss <- is.na(ord)
	if(any(miss))
	{
		message("sort found NA in ord")
		browser()
#		high.ord <- (max(ord, na.rm = TRUE) + 1):length(ord)
#		if(length(high.ord) != sum(miss))
#		{ message("high.ord error"); browser()}
#		ord0 <- ord
#		ord[miss] <- high.ord
#		if(any(is.na(ord)) || any(1:length(ord) != sort(ord)))
#		{ message("missing values of ord"); browser()}
	}
	nam <- names(x)[ord]
	object <- x[ord]
	abs1 <- abs(object@P0@Location1)
	if(verbose)
	{
		message("decreasing == ", decreasing)
		message("by.Location1 == ", by.Location1)
	}
	if(decreasing && by.Location1)
	{
		if(abs1[1] == max(abs1, na.rm = TRUE))
		{
			if(verbose) message("sorted properly by location")
		}
		else
		{ message("did not sort properly by location"); browser()}
	}
	if(any(is.na(names(object@loss0))))
	{ 
	  message("sort added NAs to names")
	  browser()
	}
	if(verbose)
		cat("checking validity of sorted ", class(object), "... ", sep = "")
	if(validObject(object))
	{
		if(verbose)
			cat("...sorted object is valid.\n")
	}
	else 
		stop("\nsort deprived loss of its validity")
	object
})

setClass("losses", representation("list", test.names = "character"))
setValidity("losses", function(object)
{
	len.ok <- length(object) == length(object@test.names)
	cla.ok <- all(sapply(object, is, class2 = "loss"))
	ok <- len.ok && cla.ok
	if(!ok) printInvalid(object)
	ok
})
setAs(from = "list", to = "losses", function(from)
{
	test.names <- if(is.character(names(from)))
		names(from)
	else
		sapply(from, annotation)
	stopifnot(is.character(test.names))
	losses(from, test.names = test.names)
})
setAs(from = "losses", to = "character", function(from)
{
	get.main <- function(x)
	{
		if(x == "t.test")
			"Student t test"
		else if(x == "wilcox.test")
			"Wilcoxon test"
		else if(x == "moderated.t.test")
			"moderated t test; linear model" # Smyth (2004)
		else if(x == "moderated.t.test.eFdr0")
			"moderated t test; theoretical null" # Efron (2004)
		else
			x
	}
	sapply(from@test.names, get.main)
})
setMethod("plot", signature(x = "losses", y = "missing"), function(x, y, call.par, ...)
{
	if(missing(call.par))
		call.par <- TRUE
	if(call.par) par(mfrow = c(2, 2))
	nam <- as(x, "character")
	for(i in 1:length(x))
		plot(x[[i]], main = nam[i], plot.predictionError = TRUE, ...)
	if(length(x) == 2)
		plot(x = x[[1]], y = x[[2]], xlab = nam[1], ylab = nam[2], call.par = FALSE, include.individual.plots = FALSE, ...)
})
setMethod("Stripchart", signature(x = "losses"), function(x, group.names, exponents, pch, circle, triangle, legend.x, text.x, strip, i0, i1, i.arbitrarize, FUN, sort.trivials.by.ratio, call.browser, verbose, foldchange.threshold, ...)
{
	if(missing(verbose)) verbose <- FALSE
	if(missing(sort.trivials.by.ratio))
	{
		sort.trivials.by.ratio <- TRUE
		message("sort.trivials.by.ratio = ", sort.trivials.by.ratio)
	}
	if(missing(FUN))
	{
		fun.name <- "identity"
		message("FUN = ", fun.name, "; consider FUN = exp_root instead.")
		FUN <- eval(parse(text = fun.name))
	}
	if(missing(strip))
	{
		strip <- TRUE
		message("strip = ", TRUE)
	}
	if(missing(exponents))
	{
		exponents <- numeric(0)
		message("exponents = ", exponents, "; consider NULL for automatic selection of exponents")
	}
	if(is.null(exponents))
	{
		min.length <- min(sapply(x, function(lo){length(lo@lossBayes)}))
		exponents <- 0:min(7, floor(log2(min.length / 100)))
	}
	stopifnot(is.numeric(exponents))
	if(length(exponents) != 0)
		ranks <- 100 * 2 ^ exponents
	get.vec <- function(lo, name, by.Location1)
	{
		message("Stripchart sorting ", class(lo))
		nsig <- sum(lo@P0 < 1, na.rm = TRUE)
		nfeatures <- max(ranks)
		if(nsig < nfeatures)
			warning(paste("Stripchart: Only ", nsig, " features have nonzero probabilities of differential expression, but ", nfeatures, " features are plotted.", sep = ""))
		if(verbose) message("\n\nget.vec calling sort for class ", class(lo))
		object <- sort(lo, by.PValue = FALSE, by.Location1 = by.Location1, verbose = verbose)
		slot(object = object, name = name)
	}
	get.elem <- function(lo, name, by.Location1)
	{	
		if(missing(by.Location1))
			by.Location1 <- FALSE #!missing(nulls.boo)
#		decreasing <- by.Location1 && !nulls.boo
		if(missing(name))
			name <- "lossBayes"
		stopifnot(is(lo, "loss"))
		elem <- if(length(exponents) == 0)
		{
			get.error <- function(name){slot(lo, name = name)}
			raw.error <- get.error(name = name)
			h0.error <- get.error(name = "loss0")
			h0.error[which(h0.error == 0)] = NaN
			absolute.error <- mean(raw.error, na.rm = TRUE) / mean(h0.error, na.rm = TRUE)
			stopifnot(length(raw.error) == length(h0.error))
			relative.error <- raw.error / h0.error
			error.mode <- hsm(relative.error, na.rm = TRUE)
			error.mean <- mean(relative.error, na.rm = TRUE)
			errors <- c(error.mode, error.mean, absolute.error)
			stopifnot(is.numeric(errors) && length(errors) == 3)
			errors
		}
		else
		{
			vec <- get.vec(lo = lo, name = name, by.Location1 = by.Location1)
			stopifnot(all(ranks <= length(vec)))
			message("Averaging error for ", name, "; ", annotation(lo), ".")
			acc <- accumulate(vec, indices = 1:max(ranks))
			stopifnot(length(acc) == max(ranks))
			acc[ranks]
		}
		FUN(elem)
	}
	which.best.loss1 <- function()
	{
		nmissing1 <- sapply(x, function(lo){sum(is.na(lo@loss1))})
		which(nmissing1 == min(nmissing1))[1]
	}
	best.i <- which.best.loss1()
	stopifnot(best.i %in% 1:length(x))
	if(missing(i0))
	{
		i0 <- best.i # 1
		message("Using ", class(x), " element i0 = ", i0, " to get predictions conditional on all null hypotheses true.")
	}
	if(missing(i1))
	{
		i1 <- best.i # i0
		message("Using ", class(x), " element i1 = ", i1, " to get predictions conditional on all null hypotheses false.")
	}
	get.trivial.elem <- function(i, name)
	{
		get.elem(x[[i]], name = name, by.Location1 = sort.trivials.by.ratio)
	}
	elem.0 <- get.trivial.elem(i = i0, name = "loss0") # get.elem(x[[i0]], name = "loss0", by.Location1 = sort.trivials.by.ratio)
	elem.1 <- get.trivial.elem(i = i1, name = "loss1")
	if(missing(i.arbitrarize))
	{
		i.arbitrarize <- numeric(0) # best.i # i0
		message("i.arbitrarize = ", i.arbitrarize)
	}
	if(missing(foldchange.threshold))
	{
		foldchange.threshold <- numeric(0) # 2 ^ c(0.5, 1, 1.5, 2)
		cat("foldchange.threshold = "); print(foldchange.threshold)
	}
	arb.lis <- if(length(foldchange.threshold) != 0)
	{
		lapply(foldchange.threshold, function(ft)
		{
			get.elem(arbitrarize(x[[i.arbitrarize]], foldchange.threshold = ft))
		})
	}
	else if(!is.null(i.arbitrarize) && is.numeric(i.arbitrarize) && length(i.arbitrarize) == 1)
	{
		arb.x <- arbitrarize(x[[i.arbitrarize]])
		elem.arb <- get.elem(arb.x)
		list(elem.arb)
	}
	else
		list()
	arb.are.ok <- function()
	{
		elem.arb.ok <- function(elem.arb)
		{
			arb.is.ok <- function(i, elem)
			{
				i1 != i.arbitrarize || all(elem.arb == elem)
			}
			ok <- is.numeric(elem.arb) && arb.is.ok(i = i0, elem = elem.0) && arb.is.ok(i = i1, elem = elem.1)
			if(!ok)
			{ message("arb.is.ok returned FALSE"); browser()}
			ok
		}
		length(arb.lis) == 0 || all(sapply(arb.lis, elem.arb.ok))
	}
	if(!arb.are.ok())
	{ message("bad arb.lis"); browser()}
	message("CALCULATING lis"); lis <- c(lapply(x, get.elem), arb.lis, list(elem.0, elem.1))
	stopifnot(all(length(lis[[1]]) == sapply(lis, length)))
	if(missing(group.names))
	{
		get.group.names <- function()
		{
			test.names <- as(x, "character")
			mle.name <- function(nulls.boo, by.Location1)
			{
				test.name <- function(i) {test.names[i]}
				prefix <- paste("all nulls ", if(nulls.boo) "true" else "false", sep = "")
				if(by.Location1)
					prefix
				else
				{
					suffix <- paste(" (sorted by ", if(nulls.boo) test.name(i0) else test.name(i1), ")", sep = "")
					paste(prefix, suffix, sep = "")
				}
			}
			arb.nam <- if(length(foldchange.threshold) != 0)
			{
				paste("fold change >", round(foldchange.threshold, 2), "(too optimistic)")
			}
			else if(length(arb.lis) == 1)
				"fold-change"
			else if(length(arb.lis) == 0)
				c()
			else
				stop("bad length(arb.lis)")
			c(test.names, arb.nam, mle.name(nulls.boo = TRUE, by.Location1 = sort.trivials.by.ratio), mle.name(nulls.boo = FALSE, by.Location1 = sort.trivials.by.ratio))
		}
		group.names <- get.group.names()
	}
	stopifnot(length(lis) == length(group.names))
	if(!missing(call.browser) && call.browser)
	{ message("Stripchart's call.browser is ", call.browser); browser()}
	if(missing(legend.x))
	{
		legend.x <- "topright"
		message("legend.x = ", legend.x)
	}
	if(missing(text.x))
	{
		text.x <- numeric(0)
		message(class(x), " text.x = ", text.x)
	}
	stopifnot(is.numeric(text.x))
	Stripchart.wrapper <- function(...){Stripchart(x = lis, group.names = group.names, text.x = text.x, ...)}
	if(length(exponents) == 0)
	{
		if(missing(pch))
			pch <- 1:length(lis[[1]])
		Stripchart.wrapper(xlab = "estimated prediction error", pch = pch, ...)
		abline(v = 1, col = "gray")
		leg <- c("relative error mode", "relative error mean", "absolute error")
		if(length(leg) != length(pch))
		{ message("pch does not match leg"); browser()}
		text.col <- gray(level = 0.7)
		legend(x = legend.x, legend = leg, pch = pch, col = text.col, text.col = text.col, bty = "n")
	}
	else if(strip)
	{
		if(missing(pch))
			pch <- as.character(exponents)
		if(missing(circle))
		{
			circle <- "4"
			message("circle = ", circle)
		}
		circle.ok <- circle %in% pch || circle %in% 1:length(pch)
		if(!circle.ok)
		{
			message("circle ", circle, " reset to default")
			circle <- character(0)
		}
		if(missing(triangle))
		{
			triangle <- "3"
			message("triangle = ", triangle)
		}
		triangle.ok <- triangle %in% pch || triangle %in% 1:length(pch)
		if(!triangle.ok)
		{
			message("triangle ", triangle, " reset to default")
			triangle <- character(0)
		}
		max.level <- 0.7
		level <- if(length(exponents) == 1)
			max.level
		else
			max.level * (max(exponents) - exponents) / max(exponents)
		if(is.na(level) || level < 0 || level > max.level)
		{ message("bad level: ", level); browser()}
		col <- gray(level = level)
		Stripchart.wrapper(xlab = "average estimated prediction error", pch = pch, col = col, circle = circle, triangle = triangle, ...)
		stopifnot(length(pch) == length(ranks))
		legend(x = legend.x, legend = paste("1", as.character(ranks), sep = "-"), pch = pch, col = col, title = "   ranks")
	}
	else
	{
		stop("non-strip Stripchart not yet implemented")
	}
})

setClass("ProbabilityEstimate", representation("Numeric", qRatio = "Scalar", lower.tail = "logical", null.hypothesis = "logical"))
setValidity("ProbabilityEstimate", function(object)
{
	is.small <- function(x){length(x) %in% c(0, 1)}
	PE.ok <- all(object <= 1.001, na.rm = TRUE) && length(object@lower.tail) <= 2 && is.small(object@null.hypothesis)
	if(!PE.ok)
	{ printInvalid(object = object); browser()}
	PE.ok
})
setMethod("annotation", signature(object = "ProbabilityEstimate"), function(object)
{
	ProbabilityEstimate.lab <- function(qRatio, lower.tail)
	{
		prefix <- "prob. of"
		suffix <- round(qRatio, 2)
		if(length(lower.tail) == 0)
			paste(prefix, "|ratio| >", suffix)
		else if(length(lower.tail) == 2)
		{
			pe.lab <- function(qRatio, lower.tail)
			{
				if(lower.tail)
					qRatio <- 1 / qRatio
				ProbabilityEstimate.lab(qRatio = qRatio, lower.tail = lower.tail)
			}
			paste(pe.lab(qRatio = qRatio, lower.tail = lower.tail[1]), "or", pe.lab(qRatio = qRatio, lower.tail = lower.tail[2]))
		}
		else if(lower.tail)
			paste(prefix, "ratio <", suffix)
		else
			paste(prefix, "ratio >", suffix)
	}
	ProbabilityEstimate.lab(qRatio = object@qRatio, lower.tail = object@lower.tail)
})
zplot <- function(x, y, ...){printGeneric("zplot"); browser()} # new 25 April 2008

#XXX|:from data.s------------------
#removeMethods("zplot")
setMethod("zplot", signature(x = "Fdr", y = "missing"), function(x, y, ...)#XXX| class"Fdr" not defined
{
	get.lab <- function(prefix){paste(prefix, "(z space)")}
	plot(zvalue(x, type = "pvalue"), zvalue(x, type = "fdr"), xlab = get.lab(prefix = "p-value"), ylab = get.lab(prefix = "local false discovery rate"), ...)
	abcol <- "gray"
	abline(v = 0, col = abcol)
	abline(h = 0, col = abcol)
})
setMethod("zplot", signature(x = "Prob0", y = "missing"), function(x, y, ...)#XXX| class"Prob0" not defined
{
	zplot(x = as(x, "Fdr"), ...)
})
#------------------
# end of file:

#Source("fdr.s")