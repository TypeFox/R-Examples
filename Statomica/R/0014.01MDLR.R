# MDL.r created by David Bickel on 14 July 2010.

#library(distr)
#library(distrDoc) # conflicted with data.s in R version 2.9.1 (2009-06-26)
#Source("data.r") # likelihood.r")
#source("distribution.r")
#source("Efron.s")

nml.type <- "normalized maximum likelihood"##XXX|: from MDL.S
#----
setClass("Family", representation(Distr.fun = "function", unknown.param.name = "character", known.param = "list", discrete = "logical", lower.unknown.param = "scalar", upper.unknown.param = "scalar"))
setValidity("Family", function(object)
{
	if(is.nothing(object))
		TRUE
	else
	{
		cla.ok <- length(object@known.param) == 0 || is.character(names(object@known.param))
		len.ok <- length(object@unknown.param.name) >= 1 && length(object@discrete) == 1
		bounds.ok <- (length(object@lower.unknown.param) == 0 && length(object@upper.unknown.param) == 0) || object@lower.unknown.param <= object@upper.unknown.param
		ok <- cla.ok && len.ok && bounds.ok
		if(!ok)
		{ printInvalid(object); browser()}
		ok
	}
})
setClass("extendedFamily", representation("Family", mle.fun = "function")) # , lower.mle = "scalar", upper.mle = "scalar"))
setValidity("extendedFamily", function(object)
{
	is.ok <- function(fun)
	{
		"x" %in% argnames(fun) # !is(fun, "likFun")
	}
	funs.ok <- is.ok(object@mle.fun)
	ok <- funs.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setClass("finiteFamily", representation("Family", unknown.param.set = "list"))
setValidity("finiteFamily", function(object)
{
	len.ok <- length(object@lower.unknown.param) == 0 && length(object@upper.unknown.param) == 0
	nam.ok <- sameSet(names(object@unknown.param.set), object@unknown.param.name)
	oks <- c(len = len.ok, nam = nam.ok)
	ok <- all(oks==T)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setClassUnion("familyContaining", c("Family"))
setMethod("family", signature(object = "Family"), function(object)
{
	object
})
setClass("Families", representation("list"))
setValidity("Families", function(object)
{
	cla.ok <- all(are(object, "Family"))
	len.ok <- length(object) > 0
	ok <- cla.ok && len.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})

setClass("error", representation("Numeric", ann = "character", estimator = "function"))
setMethod("plot", signature(x = "error", y = "error"), function(x, y, main = "", xlab, ylab, suffix = "error", add = FALSE, ...)
{
	stopifnot(length(x) == length(y))
	get.lab <- function(object){paste(annotation(object), suffix)}
	if(missing(xlab))
		xlab <- get.lab(x)
	if(missing(ylab))
		ylab <- get.lab(y)
	graph <- function(fun, ...)
	{
		mess <- function(object){message("plotting error: "); print(object)}
		mess(x)
		mess(y)
		fun(as.numeric(x), as.numeric(y), ...)
	}
	message("\n", if(add) "calling points" else "calling plot")
	if(add)
		graph(fun = points, ...)
	else
		graph(fun = plot, main = main, xlab = xlab, ylab = ylab, ...)
})
setMethod("mean", signature(x = "error"), function(x)
{
	x@.Data <- mean(as(x, "numeric"))
	x
})
setMethod("[", signature(x = "error", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@.Data <- as(x, "Numeric")[i]
  x
})

setClass("errorPair", representation(error0 = "error", error1 = "error"))
setValidity("errorPair", function(object)
{
	len.ok <- length(object@error0) == length(object@error1)
	id <- function(name){identical(slot(object@error0, name = name), slot(object@error1, name = name))}
	ann.ok <- id("ann")
	estimator.ok <- ann.ok # id("estimator")
	ok <- len.ok && ann.ok && estimator.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("hist", signature(x = "errorPair"), function(x, call.par = TRUE, call.plot = TRUE, ...)
{
	if(call.par)
		Mfrow()
	his <- function(err, xlab)
	{
		hist(err, xlab = xlab, main = annotation(err), ...)
		message(paste("\n", xlab, ":", sep = ""))
		print(stats(err))
	}
	his(x@error0, xlab = "error under null")
	his(x@error1, xlab = "error under alternative")
	if(call.plot)
		plot(x, verbose = FALSE)
})
setMethod("plot", signature(x = "errorPair", y = "missing"), function(x, y, verbose = TRUE, ...)
{
	if(verbose)
		message("consider hist(x)")
	plot(new.errorPairs(list(x)), ...)
})
coverage.root <- function(noncoverage.rate)
{
	stopifnot(is.prob(noncoverage.rate))
	1 - noncoverage.rate
}
setMethod("plot", signature(x = "numeric", y = "errorPair"), function(x, y, add = FALSE, versus.P0 = FALSE, xlab = if(versus.P0) "proportion of true nulls" else "proportion of false nulls", coverage = FALSE, root = if(coverage) coverage.root else default(TRUE, "root"), type = "l", relative = FALSE, ylab, ...)
{
	if(missing(ylab))
	{
		ylab <- if(coverage)
		{
			if(relative) "relative coverage rate?" else "coverage rate"
		}
		else
		{
			if(relative) "relative error" else "absolute error"
		}
	}
	fun <- if(add)
	{
		if(type == "l")
			lines
		else if(type == "p")
			points
		else
			stop("not my type")
	}
	else
		function(...){plot(..., xlab = xlab, ylab = ylab, type = type)}
	stopifnot(is.prob(x))
	risk <- sapply(x, function(X)
	{
		P0 <- if(versus.P0) X else 1 - X
		if(coverage)
			stopifnot(is(root, "function"))
		risk.estimate(y, P0 = P0, root = root, relative = relative)
	})
	fun(x = x, y = risk, ...)
})

setClass("errorPairs", representation("list"))
setValidity("errorPairs", function(object)
{
	ok <- all(sapply(object, is, class2 = "errorPair"))
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "numeric", y = "errorPairs"), function(x, y, versus.P0 = FALSE, legend.x, lty = 1:length(y), lwd = rep(2, length(y)), col = 1:length(y), xlim, ylim, coverage = FALSE, root = if(coverage) coverage.root else default(TRUE, "root"), relative = FALSE, call.browser = FALSE, call.par = FALSE, ann, ...)
{
	if(call.par)
		par(mfrow = c(2, 1))
	if(missing(legend.x))
	{
		legend.x <- if(relative)
			"top"
		else
		{ if(versus.P0) "topright" else "topleft"}
	}
	stopifnot(is.prob(x))
	if(missing(xlim))
		xlim <- range(x)
	if(missing(ylim))
	{
		ylim <- if(TRUE) # relative)
		{
			P0 <- if(versus.P0) max(xlim) else 1 - min(xlim)
			if(coverage)
				stopifnot(is(root, "function"))
			risks <- sapply(y, risk.estimate, root = root, P0 = P0, relative = relative)
			ymax <- max(if(call.par) 1 else 2, max(risks))
			if(call.browser)
			{
				message("risks:"); print(stats(risks))
				browser()
			}
			c(0, ymax)
		}
		else
			c(0, 1)
	}
	col <- ifelse(col == 5, "black", col)
	for(i in 1:length(y))
	{
		color <- col[i]
		message("  color: ", color)
		plot(x = x, y = y[[i]], add = i > 1, lty = lty[i], lwd = lwd[i], col = color, xlim = xlim, ylim = ylim, versus.P0 = versus.P0, relative = relative, coverage = coverage, root = root, ...)
	}
	if(!versus.P0)
		abline(v = 0.5, col = "gray")
	if(call.par) blank.plot()
	leg <- sapply(y, annotation)
	if(!missing(ann))
		leg <- ifelse(ann == "", leg, ann)
	legend(x = legend.x, legend = leg, lty = lty, lwd = lwd, col = col)
})
setMethod("plot", signature(x = "errorPairs", y = "missing"), function(x, y, length.out = 100, versus.P0 = FALSE, log = if(versus.P0) "" else "x", xlim, ...)
{
	logarithmic <- log == "x" || log == "xy"
	if(missing(xlim))
		xlim <- if(logarithmic) c(1e-4, 1) else c(0, 1)
	P0 <- if(logarithmic)
		exp(seq(log(xlim[1]), log(xlim[2]), length.out = length.out))
	else
		seq(xlim[1], xlim[2], length.out = length.out)
	stopifnot(is.prob(P0))
	plot(x = P0, y = x, xlim = xlim, log = log, versus.P0 = versus.P0, ...)
})
setMethod("[", signature(x = "errorPairs", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	lis <- as(x, "list")[i]
	nam <- names(x)[i]
  x@.Data <- lis
  names(x) <- nam
  stopifnot(validObject(x))
  x
})

setClass("posteriorP0.error", representation("error", portion.above = "Numeric", bias = "numeric", P0.bias = "scalar", ncp.bias = "scalar", score = "Numeric", null.indicator = "Numeric", unknown.param = "list", x.fun = "function"))
setValidity("posteriorP0.error", function(object)
{
	null.indicator.ok <- length(object@null.indicator) == 1 || all(object@null.indicator %in% c(0, 1)) # is.prob(object@null.indicator)
	nam.ok <- compatible(object, object@null.indicator, object@bias, object@score)
	portion.above.ok <- is.prob(object@portion.above)
	ok <- null.indicator.ok && nam.ok && portion.above.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setAs(from = "posteriorP0.error", to = "data.frame", function(from)
{
	digits <- 3
	get.vec <- function(name){round(as(slot(from, name = name), "numeric"), digits = digits)}
	data.frame(null.indicator = get.vec("null.indicator"), MSE = round(as(from, "numeric"), digits = digits), score = get.vec("score"), bias = get.vec("bias"), portion.above = get.vec("portion.above")) # , P0.bias = get.vec("P0.bias"), ncp.bias = get.vec("ncp.bias"))
})
setMethod("plot", signature(x = "posteriorP0.error", y = "posteriorP0.error"), function(x, y, suffix = "RMSE", xlim, ylim, ...)
{
	get.error <- function(object, use.null = NULL)
	{
		boo <- if(is.null(use.null))
			1:length(object)
		else if(use.null)
			1
		else
			0
		err <- sqrt(mean(as(object[object@null.indicator == boo], "error")))
		char <- if(is.null(use.null))
			""
		else
		{	if(use.null) "null" else "alt."}
		me <- try(message(annotation(object), " ", char, " RMSE: ", as(err, "numeric")))
		if(is.err(me))
		{ message("bad me"); browser()}
		err
	}
	alt.col <- "black"
	null.col <- "black"
	alt.pch <- 4
	null.pch <- 1
	if(missing(xlim) || missing(ylim))
	{
		lim <- range(get.error(x, TRUE), get.error(y, TRUE), get.error(x, FALSE), get.error(y, FALSE))
		if(missing(xlim))
			xlim <- lim
		if(missing(ylim))
			ylim <- lim
	}
	graph <- function(use.null, add, ...)
	{
		plot(x = get.error(x, use.null = use.null), y = get.error(y, use.null = use.null), add = add, col = if(use.null) null.col else alt.col, pch = if(use.null) null.pch else alt.pch, suffix = suffix, ...)
	}
	graph(use.null = TRUE, add = FALSE, xlim = xlim, ylim = ylim, ...)
	graph(use.null = FALSE, add = TRUE, ...)
	abline(a = 0, b = 1, col = "gray")
	legend(legend = c("true null hypothesis", "false null hypothesis"), x = "center", pch = c(null.pch, alt.pch), col = c(null.col, alt.col))
	print(x)
	print(y)
})
setMethod("print", signature(x = "posteriorP0.error"), function(x)
{
	digits <- 3
	datf <- as(x, "data.frame")
	nro <- 10
	min.i <- if(nrow(datf) <= nro)
		1
	else
		min(nrow(datf), max(1, max(which(x@null.indicator == 0)) - round(nro / 2)))
	max.i <- min(min.i + nro - 1, nrow(datf))
	stopifnot(min.i <= max.i)
	message("\n", annotation(x))
	print(datf[min.i:max.i, ])
	mean.stri <- function(name = "MSE", first = FALSE)
	{
		vec <- if(name == "MSE")
			as(x, "numeric")
		else
			slot(x, name = name)
		paste(if(first) "" else "; ", name, ": ", round(mean(vec), digits), sep = "")
	}
	message(mean.stri(first = TRUE), mean.stri("score"), mean.stri("bias"), mean.stri("portion.above"), "; P0.bias: ", round(x@P0.bias, digits), "; ncp.bias: ", round(x@ncp.bias, digits))
	unknown.param <- unknownParam(x)
	vec <- as(unknown.param, "numeric")
	param <- if(length(vec) == length(unknown.param))
	{
		names(vec) <- names(unknown.param)
		vec
	}
	else
		unknown.param
	print(param)
})
setMethod("mean", signature(x = "posteriorP0.error"), function(x)
{
	x@.Data <- mean(as(x, "numeric"))
	x@null.indicator <- mean(x@null.indicator)
	x@bias <- mean(x@bias)
	x@score <- mean(x@score)
	x
})
setMethod("[", signature(x = "posteriorP0.error", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@.Data <- as(x, "error")[i]
	get.vec <- function(name){slot(x, name = name)[i]}
	x@bias <- get.vec("bias")
	x@score <- get.vec("score")
	x@null.indicator <- get.vec("null.indicator")
  x
})

setClass("posteriorP0.errors", representation("list"))
setValidity("posteriorP0.errors", function(object)
{
	cla.ok <- are(object, "posteriorP0.error")
	ok <- cla.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "posteriorP0.errors", y = "missing"), function(x, y, suffix = "RMSE", xlim, ylim, col = 1:length(x), ...)
{
	stop("Consider calling ")
})
setMethod("print", signature(x = "posteriorP0.errors"), function(x)
{
	for(i in 1:length(x))
	{
		message(" \n", names(x)[i])
		print(x[[i]])
	}
})

setClass("posteriorP0.errors.lis", representation("list", nfeature = "Numeric"))
setValidity("posteriorP0.errors.lis", function(object)
{
	cla.ok <- are(object, "posteriorP0.errors")
	len.ok <- length(object) == length(object@nfeature)
	nam.ok <- all(names(object[[1]]) == sapply(object, names))
	ok <- cla.ok && len.ok && nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "posteriorP0.errors.lis", y = "missing"), function(x, y, suffix = "RMSE", ylim = numeric(0), pch, legend.i = 1, err.measure.name = c("MSE", "score", "bias", "portion.above"), ylab = err.measure.name, method.name = names(x[[1]]), legend = method.name, legend.x = "topright", ...)
{
	stopifnot(length(ylab) == length(err.measure.name))
	if(is.null(names(ylab)))
		names(ylab) <- err.measure.name
	stopifnot(all(names(ylab) == err.measure.name))
	stopifnot(length(legend) == length(method.name))
	if(missing(pch))
		pch <- 1:length(method.name)
	stopifnot(length(pch) == length(method.name))
	get.mean <- function(method, err.measure)
	{
		stopifnot(method %in% method.name)
		stopifnot(err.measure %in% err.measure.name)
		mu.hat <- sapply(1:length(x@nfeature), function(i)
		{
			errs <- x[[i]]
			assert.is(errs, "posteriorP0.errors")
			object <- errs[method][[1]]
			if(is(object, "posteriorP0.error"))
			{
				vec <- if(err.measure %in% c("MSE", "RMSE"))
					as(object, "numeric")
				else
					slot(object, name = err.measure)
				mean(vec)
			}
			else
				stop("bad get.mean")
		})
		as(if(err.measure == "RMSE") sqrt(mu.hat) else mu.hat, "numeric")
	}
	if(length(err.measure.name) > 1)
		Mfrow()
	for(err.measure in err.measure.name)
	{
		for(j in 1:length(method.name))
		{
			method <- method.name[j]
			err.vec <- get.mean(method = method, err.measure = err.measure)
			graph <- function(fun, ...){fun(x = x@nfeature, y = err.vec, pch = pch[j], ...)}
			if(j == 1)
			{
				ylim1 <- if(is.nothing(ylim))
				{
					if(err.measure == "bias") c(-1, 1) else c(0, 1)
				}
				else
					ylim
				graph(fun = plot, ylim = ylim1, xlab = "number of null hypotheses", ylab = ylab[err.measure], log = "x", ...) # "number of features"
			}
			else
				graph(fun = points)
			if(err.measure == "bias")
				abline(h = 0, col = "gray")
			else if(err.measure == "portion.above")
				abline(h = 0.5, col = "gray")
			else
				abline(h = 0, col = "gray")
		}
		if(err.measure %in% legend.i || which(err.measure == err.measure.name) %in% legend.i)
			legend(x = legend.x, legend = legend, pch = pch)
	}
})


setClass("posteriorP0.errorPair", representation(uncorrected = "posteriorP0.error", corrected = "posteriorP0.error"))
setValidity("posteriorP0.errorPair", function(object)
{
	errorPair.ok <- validObject(as(object, "errorPair"))
	ann.ok <- annotation(object@uncorrected) == "posteriorP0" && annotation(object@corrected) == "mean.posteriorP0"
	ok <- errorPair.ok && ann.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setAs(from = "posteriorP0.errorPair", to = "errorPair", function(from)
{
	get.error <- function(name, ann)
	{
		error <- as(slot(from, name = name), "error")
		if(!missing(ann))
			error@ann <- ann
		error
	}
	ann <- paste(annotation(get.error("uncorrected")), annotation(get.error("corrected")), sep = ".")
	new.errorPair(error0 = get.error("uncorrected", ann = ann), error1 = get.error("corrected", ann = ann))
})
setMethod("plot", signature(x = "posteriorP0.errorPair", y = "missing"), function(x, y, xlab = "uncorrected LFDR RMSE", ylab = "corrected LFDR RMSE", ...)
{
	plot(x = x@uncorrected, y = x@corrected, xlab = xlab, ylab = ylab, ...)
})

setMethod("sqrt", signature(x = "error"), function(x)
{
	vec <- as(x, "Numeric")
	nam <- names(x)
  x@.Data <- sqrt(vec)
  names(x) <- nam
  stopifnot(validObject(x))
  x
})
setMethod("sqrt", signature(x = "errorPair"), function(x)
{
	x@error0 <- sqrt(x@error0)
	x@error1 <- sqrt(x@error1)
  stopifnot(validObject(x))
  x
})
setMethod("sqrt", signature(x = "errorPairs"), function(x)
{
	lis <- as(x, "list")
	nam <- names(x)
  x@.Data <- lapply(lis, sqrt)
  names(x) <- nam
  stopifnot(validObject(x))
  x
})

#XXX|:changing from mle to nmle
setClass("nmle", representation("numeric", unknown.param = "list"))
setValidity("nmle", function(object)
{#XXX|:changing from mle to nmle
	ok <- if(is.nothing(object))
	{
#		if(!is.nothing(object@unknown.param))
#			warning("this nmle object is not backwards compatible")
		TRUE
	}
	else
	{
		backward.compatible <- TRUE#"ml" %in% names(object) || "mle" %in% names(object) # length(object@ml) + length(object@mle) >= 1
		nam.ok <- length(object@unknown.param) == 0 || is.character(names(object@unknown.param))
		unknown.param.ok <- !is(object@unknown.param, "unknownParams")
		backward.compatible && nam.ok && unknown.param.ok # sameNames(blank, object@unknown.param, order.sensitive = FALSE)
	}
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("P0", signature(object = "nmle"), function(object, ...)
{#XXX|:changing from mle to nmle
	object@unknown.param$P0
	
})
setAs(from = "nmle", to = "numeric", function(from)
{#XXX|:changing from mle to nmle
	vec <- from@.Data
	names(vec) <- names(from)
	vec
})
#XXX|:changing from mle to nmle
setClass("posteriorP0", representation("Numeric", estimate = "nmle", x = "numeric", family = "Family", distr0 = "univariateDistribution"))
setValidity("posteriorP0", function(object)
{
	prob.ok <- all(is.prob(object))
	len.ok <- try(length(object@x) %in% c(0, length(object)))
	if(is.err(len.ok))
		len.ok <- TRUE
	nam.ok <- sameNames(object, object@x)
	stopifnot(nam.ok)
	ok <- prob.ok && len.ok && nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setAs(from = "posteriorP0", to = "Numeric", function(from)
{
	vec <- from@.Data
	names(vec) <- names(from)
	vec
})
setAs(from = "posteriorP0", to = "numeric", function(from)
{
	as(as(from, "Numeric"), "numeric")
})
#-----------------XXX|:fails
	
setMethod("P0", signature(object = "posteriorP0"), function(object, ...)
{
	P0(object@estimate)
})
#-----------------------
#setMethod("stats", signature(object = "posteriorP0"), function(object){stats(as(object, "numeric"))})
setMethod("plot", signature(x = "posteriorP0", y = "missing"), function(x, y, xlab = "posterior probability threshold", ylab = expression(paste(proportion <= threshold)), xlim = c(0, 1), v = as.numeric(NA), ...)
{
	plot(ecdf(x), xlab = xlab, ylab = ylab, xlim = xlim, ...)
	if(is.na(v))
	{
		p0 <- try(P0(x))
		if(!is.err(p0))
			v <- p0
	}
	message("v = ", v)
	if(length(v) == 1 && !is.na(v))
		abline(v = v, col = "gray")
})
setMethod("plot", signature(x = "posteriorP0", y = "posteriorP0"), function(x, y, xlim = c(0, 1), ylim = xlim, xlab = "x posterior probability", ylab = "y posterior probability", ...)
{
	stopifnot(sameNames(x, y))
	Mfrow()
	graph <- function(post0, main)
	{
		plot(post0, main = main)
	}
	graph(x, main = xlab)
	graph(y, main = ylab)
	plot(as(x, "numeric"), as(y, "numeric"), xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
	abline(a = 0, b = 1, col = "gray")
})
setMethod("names", signature(x = "posteriorP0"), function(x){names(x@x)})

setClass("posteriorP0plus", representation("posteriorP0", estimates = "matrix"))

setClass("unknownParams", representation("list"))
setValidity("unknownParams", function(object)
{
	cla.ok <- all(sapply(object, is, class2 = "list"))
	nam.ok <- all(sapply(object, function(elem){sameNames(elem, object[[1]])}))
	ok <- cla.ok && nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})

setClass("resampled.posteriorP0", representation("matrix", unknown.params = "unknownParams", x = "matrix", original = "posteriorP0", posteriorP0.fun = "function"))
setValidity("resampled.posteriorP0", function(object)
{
	prob.ok <- is.prob(as(object, "numeric"))
	dim.ok <- all(dim(object) == dim(object@x))
	nam.ok <- all(sapply(list(rownames(object@x), names(object@original)), function(nam){length(nam) == nrow(object) && all(rownames(object) == nam)}))
	ok <- prob.ok && dim.ok && nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "resampled.posteriorP0", y = "missing"), function(x, y, file = paste("LFDR ", Sys.Date(), " ", round(runif(1) * 1000), ".pdf", sep = ""), ...)
{
	point0 <- x@original
	expected0 <- as(x, "posteriorP0")
	plot(x = point0, y = expected0, xlab = "original LFDR point estimate", ylab = "confidence mean LFDR", ...)
	pdf(file = file)
	print(file)
	for(i in 1:nrow(x))
	{
		hist(x[i, ], xlab = "LFDR estimate", main = rownames(x)[i])
		abline(v = x@original[i], col = "gray")
	}
	dev.off()
})
setMethod("mean", signature(x = "resampled.posteriorP0"), function(x, method = "p", ptrans = pnorm, qtrans = qnorm)
{
	stopifnot(length(method) == 1)
	post0 <- x@original
	get.mean <- function(bootstrap.estimate, original.estimate)
	{
		expected <- if(method == "p") # percentile
			mean(bootstrap.estimate) # original.estimate, ptrans, and qtrans ignored
		else if(method == "t") # percentile t
		{
			stop("needs variance estimator, perhaps from double bootstrapping (shao:tu:1995, p. 131)")
		}
		else if(method == "bcp") # bias-corrected percentile
		{
			if(!(length(bootstrap.estimate) > 1 && length(original.estimate) == 1))
			{ message("bad length for mean of resampled.posteriorP0 object"); browser()}
			params <- BC.percentile(xstar = bootstrap.estimate, point = original.estimate, ptrans = ptrans, qtrans = qtrans)
			mean(params)
		}
		else if(method == "bc") # bias-corrected
		{
			bootstrap.q <- mean(qtrans(bootstrap.estimate))
			original.q <- qtrans(original.estimate)
			bias.q <- bootstrap.q - original.q
			stopifnot(length(bias.q) == 1)
			if(!is.finite(bias.q))
				bias.q <- 0
			ptrans(original.q - bias.q) # == 2 * original.q - bootstrap.q
		}
		else
			stop("method not recognized")
		if(!is.prob(expected))
		{ message("bad expected: "); print(stats(expected)); browser()}
		expected
	}
	expected0 <- sapply(1:nrow(x), function(i){get.mean(bootstrap.estimate = x[i, ], original.estimate = post0[i])})
	stopifnot(length(expected0) == length(post0))
	names(expected0) <- names(post0)
	post0@.Data <- Numeric(expected0)
	if(length(x@unknown.params) == 1)
		post0@estimate@unknown.param <- x@unknown.params[[1]]
	else
		warning("retaining unknown.param of x@original")
	names(post0) <- names(expected0)
	post0
})
setAs(from = "resampled.posteriorP0", to = "posteriorP0", function(from)
{
	post0 <- mean(from)
	assert.is(post0, "posteriorP0")
	post0
})


setClass("Weight", representation(focus.weight = "scalar", incidental.weight = "numeric", lower.bound = "scalar", upper.bound = "scalar"))
setValidity("Weight", function(object)
{
	if(is.nothing(object))
		TRUE
	else
	{
		len.ok <- length(object@focus.weight) == 1 && length(object@lower.bound) == 1 && length(object@upper.bound) == 1
		order.ok <- object@lower.bound <= object@upper.bound
		w <- sum(c(object@focus.weight, object@incidental.weight))
		in.bounds <- w >= object@lower.bound && w <= object@upper.bound
		focused <- all(object@focus.weight >= object@incidental.weight)
		oks <- c(len = len.ok, nam = order.ok, bounds = in.bounds) #, focused = focused)
		ok <- all(oks==T)
		if(!ok)
		{ printInvalid(object); browser()}
		ok
	}
})
setClass("Weights", representation("list"))
setValidity("Weights", function(object)
{
	cla.ok <- all(are(object, "Weight"))
	len.ok <- all(ncomparison(object) == sapply(object, ncomparison))
	ok <- cla.ok && len.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})

setClass("likFun", representation("function", Distr.fun = "function", unknown.param.name = "character", known.param = "list", x = "numeric", base = "Scalar", weight = "Weight", arglis = "list"))
setValidity("likFun", function(object)
{
	cla.ok <- length(object@known.param) == 0 || is.character(names(object@known.param))
	len.ok <- length(object@base) == 1 && length(object@x) >= 1 && length(object@unknown.param.name) >= 1
	fun.ok <- !is(object@Distr.fun, "likFun") && !("x" %in% argnames(object@Distr.fun))
	arglis.ok <- length(object@arglis) == 0 || is.character(names(object@arglis))
	ok <- cla.ok && len.ok && object@base > 0 && fun.ok && arglis.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setAs(from = "likFun", to = "Weight", function(from)
{
	from@weight
})
setMethod("plot", signature(x = "likFun", y = "missing"), function(x, y, xlim, length.out = 1e2, log = "", xlab, x.converges, ...)
{
	ylab <- if(x@base == 2)
		"bits for alternative hypothesis"
	else
		paste("likelihood (base ", x@base, ")", sep = "")
	if(missing(xlim))
	{
		xlim <- if(length(x@x) == 1)
		{
			div <- x@x / 2
			mult <- 2 * x@x
			lower <- min(div, mult)
			upper <- max(div, mult)
			stopifnot(lower < upper)
			c(lower, upper)
		}
		else
			stop("cannot assign xlim")
	}
	param.value <- seq(from = xlim[1], to = xlim[2], length.out = length.out)
	if(missing(xlab))
		xlab <- x@unknown.param.name
	fun <- as(x, "function")
	y <- if(x.converges)
		function(...){fun(..., x.converges = x.converges)}
	else
		fun
	plot(x = param.value, y = y, xlim = xlim, xlab = xlab, ylab = ylab, log = log, ...)
	abline(h = 0, col = "gray")
	if(is.finite(x@x))
		abline(v = x@x, col = "gray")
})

setClass("lik.ratio", representation("numeric", type = "character", base = "Scalar", unknown.param.set = "list", estimate = "numeric"))
setValidity("lik.ratio", function(object)
{
	len.ok <- length(object@type) == 1 && length(object@estimate) %in% c(0, 1, length(object))
	base.ok <- if(length(object@base) == 0)
		all(object >= 0, na.rm = TRUE)
	else if(length(object@base) == 1)
		object@base > 0
	ok <- len.ok && base.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
lik.ratio.type <- function(object,nml.type=nml.type)
{#XXX:added ,nml.type as input 2014
	assert.is(object, "lik.ratio")
	if(object@type == nml.type) "NML" else object@type
}
lik.ratio.lab <- function(object, hindsight = numeric(0), long.lab = FALSE, favoring.null = TRUE, base)
{
	assert.is(object, "lik.ratio")
	use.regret <- !is.nothing(hindsight)
	clab <- if(favoring.null) "favoring null" else "favoring alternative"
	lab <- paste(lik.ratio.type(object), if(use.regret) "regret of alternative" else clab)
	if(long.lab)
	{
		lab <- paste(lab, " (", if(base == 2) "bits" else paste("log", base, sep = ""), ")", sep = "")
		if(length(object@unknown.param.set) > 0)
		{
#			lab <- expression(paste(lab, " |", Omega, "| = ", length(object@unknown.param.set), sep = "")) # expression(paste("Phase Angle ", phi))
			lab <- paste(lab, " |param. set| = ", length(object@unknown.param.set), sep = "")
			print(lab)
		}
	}
	lab
}
setMethod("plot", signature(x = "lik.ratio", y = "missing"), function(x, y, xlab = "estimate", ...)
{
	X <- x@estimate
	if(length(x) == length(X))
		plot(x = X, y = x, xlab = xlab, ...)
	else
		plot(as(x, "numeric"), ...)
})
setMethod("plot", signature(x = "numeric", y = "lik.ratio"), function(x, y, hindsight = numeric(0), base = 10, long.lab = TRUE, xlab, ylab, main, favoring.null = TRUE, add = FALSE, ...)
{
	get.cl <- function(object){codelength(object = object, hindsight = hindsight, base = base, favoring.null = favoring.null)}
	if(missing(xlab))
		xlab <- "statistic"
	if(missing(ylab))
		ylab <- lik.ratio.lab(object = y, hindsight = hindsight, long.lab = long.lab, favoring.null = favoring.null, base = base) # paste(, " (base ", base, ")", sep = "")
	if(missing(main))
		main <- lik.ratio.type(y)
	ycl <- get.cl(y)
	plot(x = x, y = ycl, xlab = xlab, ylab = ylab, main = main, ...)
})
setMethod("plot", signature(x = "lik.ratio", y = "lik.ratio"), function(x, y, hindsight = numeric(0), base = 10, long.lab = TRUE, xlab, ylab, favoring.null = TRUE, add = FALSE, ...)
{
	get.cl <- function(object){codelength(object = object, hindsight = hindsight, base = base, favoring.null = favoring.null)}
	xtype <- lik.ratio.type(x)
	ytype <- lik.ratio.type(y)
	xcl <- get.cl(x)
	ycl <- get.cl(y)
	stopifnot(sameNames(xcl, ycl))
	get.lab <- function(object)
	{
		lik.ratio.lab(object = object, hindsight = hindsight, long.lab = long.lab, favoring.null = favoring.null, base = base)
	}
	graph <- function(fun, ...)
	{
		fun(x = xcl, y = ycl, ...)
	}
	if(add)
		graph(fun = points, ...)
	else
		graph(fun = plot, xlab = if(missing(xlab)) get.lab(x) else xlab, ylab = if(missing(ylab)) get.lab(y) else ylab, ...)
	abline(a = 0, b = 1, col = "gray")
# plot(log2(Hml10sho)-log2(Nml10sho), log2(Hml10sho)-log2(Bf10sho), xlab = "NML regret", ylab = "BDD Bayes factor regret");
})

setClass("lik.ratios", representation("list"))
setValidity("lik.ratios", function(object)
{
	len.ok <- length(object) >= 1
	cla.ok <- are(object, "lik.ratio")
	base.ok <- all(sapply(object, base) == base(object))
	nam.ok <- all(sapply(object, function(elem){sameNames(elem, object[[1]])}))
	oks <- c(len.ok = len.ok, cla.ok = cla.ok, base.ok = base.ok, nam.ok = nam.ok)
	ok <- all(oks==T)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("names", signature(x = "lik.ratios"), function(x)
{
	names(x[[1]])
})
setMethod("plot", signature(x = "xprnSetObject", y = "lik.ratios"), function(x, y, hindsight = numeric(0), base = 2, long.lab = FALSE, legend.x = "topleft", multiple, difference, FUN = icv, legend.legend, favoring.null = TRUE, ylab, abline.args, ...)
{
	stopifnot(sameNames(x, y))
	cls <- codelength(object = y, hindsight = hindsight, base = base, favoring.null = favoring.null) # changed 111022
	use.regret <- !is.nothing(hindsight)
	xlab <- if(identical(FUN, icv))
		"sample inverse coefficient of variation"
	else if(identical(FUN, mean))
		"expression ratio estimate"
	else
		"estimate"
	suffix <- if(base == 2)
		""
	else
		paste("(", if(base == 2) "bits" else paste("base", base), ")", sep = "")
	clab <- if(favoring.null) "codelength favoring null" else "information favoring alternative"
	if(missing(ylab))
		ylab <- paste(if(use.regret) "regret of alternative" else clab, suffix)
	print(c(favoring.null, ylab))
	if(missing(multiple))
		multiple <- length(y) <= 4
	ratio.estimate <- if(identical(FUN, mean))
		exp(stat(x, FUN = FUN))
	else
		stat(x, FUN = FUN)
	log <- if(identical(FUN, mean))
		"x"
	else
		""
	if(missing(abline.args))
	{
		abline.args <- list(h = 0, col = "gray")
		if(identical(FUN, mean))
			abline.args <- c(list(v = 1), abline.args)
		if(identical(FUN, icv))
			abline.args <- c(list(v = 0), abline.args)
	}
	if(missing(legend.legend))
		legend.legend <- type(y)
	plot(x = ratio.estimate, y = cls, legend.legend = legend.legend, legend.x = legend.x, xlab = xlab, ylab = ylab, multiple = multiple, abline.args = abline.args, log = log, ...) # sapply(y, slot, name = "type")
	if(missing(difference))
		difference <- default(length(y) == 2, "difference")
	if(difference)
	{
		plot(x = ratio.estimate, y = cls[[1]] - cls[[2]], xlab = xlab, ylab = paste("difference in", ylab), log = log)
		do.call(abline, abline.args)
	}
})

setClass("statistic", representation("numeric"))

setClass("min.redundancy", representation("scalar", lr = "lik.ratio", family = "Family", unknown.param.1 = "list", unknown.param.2 = "list"))
setValidity("min.redundancy", function(object)
{
	len.ok <- length(object@lr) >= 1
	nam.ok <- sameNames(object@unknown.param.1, object@unknown.param.2) && sameSet(object@family@unknown.param.name, names(object@unknown.param.2))
	ok <- len.ok && nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
##XXX|: from MDL.S
error <- function(object, P0)
{
	assert.is(object, "errorPair")
	P0 <- Scalar(P0)
	stopifnot(is.prob(P0))
	err <- object@error0
	err@.Data <- P0 * object@error0 + (1 - P0) * object@error1
	stopifnot(validObject(err))
	assert.is(err, "error")
	err
}

new.errorPair <- function(error0, error1)
{
	new("errorPair", error0 = error0, error1 = error1)
}
new.errorPairs <- function(object)
{
	new("errorPairs", object)
}
# near end of file:

#source(file = "MDL.s") # functions

