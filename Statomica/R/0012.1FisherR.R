# Fisher.r created by David Bickel on 23 April 2011.

# library(multtest) # not needed since you have independent.Sidak

#source("data.r")
#source("congruity.r")

setClass("bias.corrected.zvalue", representation("numeric", uncorrected = "numeric"))
setValidity("bias.corrected.zvalue", function(object)
{
	nam.ok <- sameNames(object, object@uncorrected)
	ok <- nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "bias.corrected.zvalue", y = "missing"), function(x, y, xlab = "z", ylab = "corrected z", ...)
{
	ord <- order(x@uncorrected)
	prepare.z <- function(z)
	{
		as(z, "numeric")[ord]
	}
	plot(x = prepare.z(x@uncorrected), y = prepare.z(x), xlab = xlab, ylab = ylab, ...)
	abline(a = 0, b = 1, col = "gray")
})
setMethod("[", signature(x = "bias.corrected.zvalue", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@.Data <- as(x, "numeric")[i]
	x@uncorrected <- x@uncorrected[i]
  x
})

setClass("bias.corrected.pvalue", representation("numeric", uncorrected = "Numeric", ranks = "Numeric")) # "Numeric" caused names problem
setValidity("bias.corrected.pvalue", function(object)
{
	prob.ok <- 1#is.prob(object@.Data) && is.prob(object@uncorrected)
	nam.ok <- sameNames(object, object@uncorrected)
	ranks.ok <- is.nothing(object@ranks) || (sameNames(object, object@ranks) && all(1 <= object@ranks & object@ranks <= length(object)))
	oks <- c(prob.ok = prob.ok, nam.ok = nam.ok, ranks.ok = ranks.ok)
	ok <- all(oks==T)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "bias.corrected.pvalue", y = "missing"), function(x, y, xlim, ylim, add = FALSE, log, xlab = "p-value", ylab = "corrected p-value", call.browser = FALSE, ...)
{
	if(call.browser)
	{ message("call.browser == TRUE"); browser()}
	pval <- try(x@uncorrected)
	if(is.err(pval))
	{ message("bad pval"); browser()}
	if(missing(xlim))
		xlim <- range(c(x, pval))
	if(missing(ylim))
		ylim <- xlim
	if(missing(log))
		log <- if(all(c(xlim, ylim) > 0)) "xy" else ""
	ord <- order(pval)
	prepare.p <- function(p)
	{
		as(p, "numeric")[ord]
	}
	graph <- function(fun, ...)
	{
		fun(x = prepare.p(pval), y = prepare.p(x), ...)
	}
	if(add)
		graph(fun = points, ...)
	else
		graph(fun = plot, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, log = log, ...)
	abline(a = 0, b = 1, col = "gray")
})
setMethod("plot", signature(x = "bias.corrected.pvalue", y = "numeric"), function(x, y, xlim = c(0, 1), ylim = c(0, 1), xlab = "selection-corrected p", ylab = "achieved FWER", ...)
{
	stopifnot(length(x) == length(y))
	plot(x = as(x, "numeric"), y = y, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
	abline(a = 0, b = 1, col = "gray")
})
setMethod("[", signature(x = "bias.corrected.pvalue", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@.Data <- as(x, "numeric")[i]
	x@uncorrected <- x@uncorrected[i]
  x
})
setAs(from = "bias.corrected.pvalue", to = "numeric", function(from)
{
	vec <- from@.Data
	names(vec) <- names(from)
	vec
})

setClass("bias.corrected.certainty", representation(lower = "Numeric", upper = "Numeric", nonmonotonic = "bias.corrected.pvalue"))
setValidity("bias.corrected.certainty", function(object)
{
	prob.ok <- is.prob(object@lower) && is.prob(object@upper)
	len.ok <- all(length(object) == c(length(object@lower), length(object@upper)))
	nam.ok <- sameNames(object@lower, object@upper, object@nonmonotonic)
	lower.ok <- all(object@lower <= object@nonmonotonic)
	upper.ok <- all(object@nonmonotonic <= object@upper)
	oks <- c(prob.ok, len.ok, nam.ok, lower.ok, upper.ok)
	ok <- all(oks==T)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("[", signature(x = "bias.corrected.certainty", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@lower <- x@lower[i]
	x@upper <- x@upper[i]
	x@nonmonotonic <- x@nonmonotonic[i]
  x
})
setMethod("length", signature(x = "bias.corrected.certainty"), function(x)
{
	length(x@nonmonotonic)
})
setMethod("Rank", signature(object = "bias.corrected.certainty"), function(object, ...)
{
	nonmon <- object@nonmonotonic
	if(is.nothing(list(...)) && !is.nothing(nonmon@ranks))
		nonmon@ranks
	else
		rank(nonmon@uncorrected, ...)
})
setMethod("plot", signature(x = "bias.corrected.certainty", y = "missing"), function(x, y, xlim, ylim, pch = 3, xlab = "p-value", ylab = "adjusted p-value", fwer.pch = 4, fdr.pch = 1, nfeature = length(x), call.browser = FALSE, ...)
{
	message("\n\nplotting ", class(x), " on ", date())
	if(!is.nothing(fwer.pch))
		achieved.fwer <- independent.Sidak(x@nonmonotonic@uncorrected)
	if(!is.nothing(fdr.pch))
		achieved.FDR <- estimated.FDR(x@nonmonotonic@uncorrected)
	i <- if(nfeature == length(x))
		1:nfeature
	else
		Rank(x) <= nfeature
	x <- x[i]
	if(missing(xlim))
	{
		xlim <- if(missing(nfeature)) c(0, 1) else range(x@nonmonotonic@uncorrected)
	}
	if(missing(ylim))
	{
		ylim <- if(missing(nfeature)) c(0, 1) else c(0, max(c(x@lower, x@upper)))
	}
	ord <- order(x@nonmonotonic@uncorrected)
	prepare.p <- function(p)
	{
		as(p, "numeric")[ord]
	}
	graph <- function(fun, y, ...)
	{
		fun(x = prepare.p(x@nonmonotonic@uncorrected), y = prepare.p(y), ...)
	}
	first.graph <- function(y, col, pch = 1)
	{
		graph(fun = plot, y = y, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, col = col, pch = pch, ...)
		abline(a = 0, b = 1, col = "gray")
	}
	if(is.nothing(pch))
	{
		first.graph(y = x@nonmonotonic, col = "gray")
		graph(fun = lines, y = x@lower, col = "gray")
		graph(fun = lines, y = x@upper, col = "gray")
		graph(fun = lines, y = as(x, "bias.corrected.pvalue"), col = "black")
	}
	else
		first.graph(y = as(x, "bias.corrected.pvalue"), col = "black", pch = pch)
	if(!is.nothing(fwer.pch))
	{
		achieved.fwer <- achieved.fwer[i]
		message("\nmean achieved FWER: ", round(mean(achieved.fwer), 3))
		graph(fun = points, y = achieved.fwer, pch = fwer.pch)
	}
	if(!is.nothing(fdr.pch))
	{
		achieved.FDR <- achieved.FDR[i]
		message("\nmean achieved FDR: ", round(mean(achieved.FDR), 3))
		graph(fun = points, y = achieved.FDR, pch = fdr.pch)
	}
	if(!is.nothing(pch) && !is.nothing(fwer.pch) && !is.nothing(fdr.pch))
		legend(x = "bottomright", legend = c("evidential p", "achieved FWER", "achieved FDR"), pch = c(pch, fwer.pch, fdr.pch))
	if(!is.nothing(fdr.pch) && call.browser)
		plot(x@nonmonotonic@uncorrected, achieved.FDR)
	message("\n  finished plotting ", class(x), " on ", date())
	if(call.browser)
		browser()
}) # end "plot", signature(x = "bias.corrected.certainty", y = "missing")
setAs(from = "bias.corrected.certainty", to = "bias.corrected.pvalue", function(from)
{
	bcp <- from@nonmonotonic
	bcp@.Data <- pmean(from@lower, from@upper)
	names(bcp) <- names(from)
	bcp
})
setAs(from = "bias.corrected.certainty", to = "numeric", function(from)
{
	as(as(from, "bias.corrected.pvalue"), "numeric")
})


# near end of file:

#source(file = "Fisher.s") # functions

