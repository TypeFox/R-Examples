# cdf.r (formerly fiducial.r) created by David Bickel on 20 April 2009.

#source("data.r") # See Jeffreys.r, Bayes.r.

setClass("CDF", representation("function", min_param = "scalar", max_param = "scalar", param.name =  "character", type = "character")) # Cf. "pposterior" of logical.r
setValidity("CDF", function(object)
{
	ok <- object@min_param <= object@max_param && length(object@param.name) == 1 && length(object@type) >= 1
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
#setMethod("plot", signature(x = "CDF", y = "scalar"), function(x, y, xlab = y@param.name, ylab, reduced = default(FALSE, "CDF reduced"), ...) #
#{
#
#})
setMethod("plot", signature(x = "numeric", y = "CDF"), function(x, y, xlab = y@param.name, ylab, reduced = default(FALSE, "CDF reduced"), ...) #
{
	type <- y@type[1]
	if(sum(names(list(...)) %in% "ylim") > 1) { message("ylim in ... (", class(x), ", ", class(y), ") more than once"); browser()}
	if(!reduced) # call.par = TRUE,
	{
		par(mfrow = c(2,2))
		plot(x = x, y = as(deriv(y), "function"), xlab = xlab, ylab = paste(type, "density"))
	}
	graph <- function(fun, ylab)
	{
		plot(x = x, y = sapply(x, fun), xlab = xlab, ylab = ylab, ...) # sapply(x, fun) reinstated 19 December 2010?
	}
	if(missing(ylab))
		ylab <- paste("cumulative", type)
	graph(fun = as(y, "function"), ylab = ylab)
})
plot.CDF <- function(x, y, xlim = stop("you have to specify xlim"), ...)
{
	message('Calling "plot", signature(x = "numeric", y = "CDF")')
	plot(x = seq(xlim[1], xlim[2], length.out = 50), y = x, reduced = TRUE, type = "l", ...)
}
setMethod("plot", signature(x = "CDF", y = "missing"), function(x, y, ...)
{
	stop('try Plot or "plot", signature(x = "numeric", y = "CDF")')
})
setMethod("Plot", signature(x = "CDF", y = "missing"), plot.CDF)
setIs("CDF", "Function") # for "plot", signature(x = "Function", y = "missing") of data.r, but does not work as of 19 December 2010

setClass("CDFs", representation("list", x = "numeric", size = "numeric")) # Cf. "pposterior" of posterior.s
setValidity("CDFs", function(object)
{
	cla.ok <- all(sapply(object, is, class2 = "CDF"))
	len.ok <- length(object) == length(object@size)
	is.ok <- function(name)
	{
		vec <- sapply(object, slot, name = name)
		all(vec[1] == vec)
	}
	slo.oks <- c(min_param = is.ok("min_param"), max_param = is.ok("max_param"), param.name = is.ok("param.name"))
	slo.ok <- all(ifelse(is.na(slo.oks), FALSE, slo.oks))
	ok <- cla.ok && len.ok && slo.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "CDFs", y = "missing"), function(x, y, file = paste("CDFs", Sys.Date(), max(x@size), "pdf", sep = "."), reduced = default(FALSE, "CDFs reduced"), ...)
{
	pdf(file = file)
	if(reduced)
		par(mfrow = c(2, 2))
	for(i in 1:length(x))
	{
		cdf <- x[[i]]
		size <- x@size[i]
		main <- paste(paste(cdf@type, collapse = "; "), size)
		message("\n", i, " ", main, "; file: ", file)
		plot(cdf, main = main, sub = file, reduced = reduced, ...)
	}
	dev.off()
})
setAs(from = "CDFs", to = "CDF", function(from)
{
	cdf <- Mean(from)
	assert.is(cdf, "CDF")
	cdf
})
setAs(from = "CDFs", to = "functions", function(from)
{
	lis <- lapply(from, as, Class = "function")
	functions(lis)
})

# objects moved from cdf.r to interval.r on 24 December 2010.

# files from logical.r moved to posterior.r on 7 July 2010.


# near end of file:

#source(file = "cdf.s") # functions











