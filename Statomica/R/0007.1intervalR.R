# interval.r created by David Bickel on 24 December 2010.

#Source("congruity.r")
#Median <- function(x, ...){} # generic function for methods#Marta Oct 2013
setMethod("Median", signature(x = "numeric"), function(x,...){median(x = x, ...)})
setClass("posteriorInterval", representation(lower = "numeric", upper = "numeric"))
setValidity("posteriorInterval", function(object)
{
	ok <- if(is.nothing(object))
		TRUE
	else
	{
		len.ok <- length(object@lower) == length(object@upper)
		nam.ok <- sameNames(length(object@lower), length(object@upper))
		na.ok <- !any(is.na(c(object@lower, object@upper)))
		ord.ok <- na.ok && all(object@lower <= object@upper)
		oks <- c(len = len.ok, nam = nam.ok, na = na.ok, ord = ord.ok)
		all(oks)
	}
	if(is.na(ok) || !ok)
	{ printInvalid(object); print(oks); browser()}
	ok
})
setMethod("length", signature(x = "posteriorInterval"), function(x)
{
	length(x@lower)
})
setMethod("names", signature(x = "posteriorInterval"), function(x)
{
	names(x@lower)
})
setMethod("Combine", signature(x = "posteriorInterval", y = "posteriorInterval"), function(x, y, ...)
{
	ci <- x
	ci@lower <- Combine(x@lower, y@lower)
	ci@upper <- Combine(x@upper, y@upper)
	ci
})
setMethod("[", signature(x = "posteriorInterval", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@lower <- x@lower[i]
	x@upper <- x@upper[i]
  x
})
setAs(from = "list", to = "posteriorInterval", function(from)
{
	if(is(from, "inverseCDF") || is(from, "inverseCDFs"))
		posteriorInterval(object = from)
	else if(all(sapply(from, is, class2 = "posteriorInterval")))
	{
		mat <- sapply(from, function(int)
		{
			assert.is(int, "posteriorInterval")
			stopifnot(length(int) == 1)
			c(int@lower, int@upper)
		})
		stopifnot(nrow(mat) == 2)
		colnames(mat) <- names(from)
		as(mat, "posteriorInterval")
	}
	else
	{ message("cannot yield posteriorInterval"); browser()}
})
setAs(from = "matrix", to = "posteriorInterval", function(from)
{
	stopifnot(nrow(from) == 2)
	lower <- from[1, ]
	upper <- from[2, ]
	names(lower) <- colnames(from)
	new.posteriorInterval(lower = lower, upper = upper)
})
setMethod("plot", signature(x = "posteriorInterval", y = "posteriorInterval"), function(x, y, base = default(2, "base"), xlab = "x width", ylab = "y width", ...)
{
	w <- function(obj){width(obj, base = base)}
	wx <- w(x)
	wy <- w(y)
	lim <- range(c(wx, wy))
	add.suffix <- function(lab)
	{
		suffix <- paste("(log", base, ")", sep = "")
		paste(lab, suffix)
	}
	if(!is.nothing(base))
	{
		xlab <- add.suffix(xlab)
		ylab <- add.suffix(ylab)
	}
	plot(wx, wy, xlab = xlab, ylab = ylab, ...) # xlim = lim, ylim = lim,
	abline(a = 0, b = 1, col = "gray")
})

setClass("posteriorIntervalPlus", representation("posteriorInterval", alternative = "logical", param = "numeric", p1 = "Scalar", p2 = "Scalar"))
setValidity("posteriorIntervalPlus", function(object)
{
	ok <- if(is.nothing(object))
		FALSE # change this to TRUE if you want to allow blank posteriorIntervalPlus
	else
	{
		nam.ok <- sameNames(object, object@alternative, object@param) #, object@p1, object@p2)
		oks <- c(p1 = is.prob(object@p1), p2 = is.prob(object@p2), nam = nam.ok)
		all(oks)
	}
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "posteriorIntervalPlus", y = "posteriorIntervalPlus"), function(x, y, base = default(numeric(0), "base"), ...)
{
	plot(x = as(x, "posteriorInterval"), y = as(y, "posteriorInterval"), base = base, ...)
})
setMethod("Plot", signature(x = "posteriorIntervalPlus", y = "posteriorIntervalPlus"), function(x, y, xlab = "x error", ylab = "y error", ...)
{
	Mfrow()
	get.err <- function(object)
	{
		get.location <- function(name){slot(object, name = name)}
		location <- get.location("lower")
		stopifnot(all(location == get.location("upper")))
		location - object@param
	}
	plot(get.err(x), get.err(y), xlab = xlab, ylab = ylab, ...)
})
setMethod("hist", signature(x = "posteriorIntervalPlus"), function(x, y, exclusion = 0, call.par = FALSE, ...)
{
	if(call.par)
	{
		Mfrow()
		his <- function(theta)
		{
			get.point <- function(name){slot(x, name = name)[x@param %in% theta]} # new.log( , base = 2
			point <- get.point("lower")
			stopifnot(all(point == get.point("upper")))
			hist(x = point, xlab = "estimated mean (log2)", ylab = "number of simulations", main = paste("mean:", signif(theta, digits = 3), ...))
			abline(v = theta, col = "gray")
		}
		all.theta <- unique(x@param)
		for(theta in setdiff(all.theta, exclusion))
			his(theta)
		his(all.theta)
	}
	get.location <- function(name){slot(x, name = name)}
	location <- get.location("lower")
	stopifnot(all(location == get.location("upper")))
	hist(location - x@param, xlab = "estimate error", ylab = "number of simulations", ...)
})
setMethod("Combine", signature(x = "posteriorIntervalPlus", y = "posteriorIntervalPlus"), function(x, y, ...)
{
	stopifnot(x@p1 == y@p1)
	stopifnot(x@p2 == y@p2)
	ci <- try(Combine(as(x, "posteriorInterval"), as(y, "posteriorInterval")))
	if(is.err(ci))
	{ message("bad ci 1"); browser()}
	alternative <- Combine(x@alternative, y@alternative)
	param <- Combine(x@param, y@param)
	new.posteriorIntervalPlus(ci, alternative = alternative, param = param, p1 = x@p1, p2 = x@p2)
})
setMethod("[", signature(x = "posteriorIntervalPlus", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	ci <- try(as(x, "posteriorInterval")[i])
	if(is.err(ci))
	{ message("bad ci 2"); browser()}
	alternative <- x@alternative[i]
	param <- x@param[i]
	new.posteriorIntervalPlus(ci, alternative = alternative, param = param, p1 = x@p1, p2 = x@p2)
})


setClass("confidencePerformance", representation(coverage = "Numeric", coverage0 = "Numeric", coverage1 = "Numeric"))
setValidity("confidencePerformance", function(object)
{
	nam.ok <- sameNames(object@coverage, object@coverage0, object@coverage1)
	prob.ok <- all(sapply(list(object@coverage, object@coverage0, object@coverage1), is.prob))
	oks <- c(prob = prob.ok, nam = nam.ok)
	ok <- all(oks)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})


setClass("inverseCDF", representation("function", min_param = "scalar", max_param = "scalar", param.name =  "character", type = "character")) # Cf. "pposterior" of logical.r
setValidity("inverseCDF", function(object)
{
	ok <- object@min_param <= object@max_param && length(object@param.name) == 1 && length(object@type) >= 1
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "numeric", y = "inverseCDF"), function(x, y, xlab, ylab = y@param.name, call.browser = FALSE, plot.CDF = TRUE, ...)
{
	if(call.browser)
	{ message("calling the browser, as requested"); browser()}
	type <- y@type[1]
	if(sum(names(list(...)) %in% "ylim") > 1) { message("ylim in ... (", class(x), ", ", class(y), ") more than once"); browser()}
	graph <- function(fun, xlab)
	{
		y <- sapply(x, fun)
		if(!any(is.finite(y)))
		{ message("bad inversion of the p-value function"); browser()}
		if(plot.CDF)
			plot(y = x, x = y, ylab = xlab, xlab = ylab, type = "l", ...)
		else
			plot(x = x, y = y, xlab = xlab, ylab = ylab, type = "l", ...) # fun
	}
	if(missing(xlab))
		xlab <- paste("cumulative", type)
	graph(fun = as(y, "function"), xlab = xlab)
})
plot.inverseCDF <- function(x, y, length.out = 50, ...)
{
	message('Calling "plot", signature(x = "numeric", y = "inverseCDF")')
	from <- 1 / (2 * length.out)
	plot(x = seq(from, 1 - from, length.out = length.out), y = x, ...)
}
setMethod("plot", signature(x = "inverseCDF", y = "missing"), plot.inverseCDF)
setMethod("Plot", signature(x = "inverseCDF", y = "missing"), plot.inverseCDF)

setClass("inverseCDFs", representation("list", size = "numeric", mass0 = "Numeric"))
setValidity("inverseCDFs", function(object)
{
	cla.ok <- all(sapply(object, is, class2 = "inverseCDF"))
	len.ok <- is.nothing(object@size) || length(object) == length(object@size)
	mass0 <- object@mass0
	mass0.ok <- is.nothing(mass0) || (sameNames(mass0, object) && all(mass0 <= 1.0001))
	is.ok <- function(name)
	{
		vec <- sapply(object, slot, name = name)
		vec <- vec[vec != "no parameter"]
		is.nothing(vec) || all(vec[1] == vec)
	}
	slo.oks <- c(min_param = is.ok("min_param"), max_param = is.ok("max_param"), param.name = is.ok("param.name"))
	slo.ok <- all(ifelse(is.na(slo.oks), FALSE, slo.oks))
	oks <- c(cla = cla.ok, len = len.ok, slo = slo.ok, mass0 = mass0.ok)
	ok <- all(oks)
	if(!ok)
	{ printInvalid(object); print(oks); browser()}
	ok
})
setAs(from = "inverseCDFs", to = "list", function(from)
{
	lis <- from@.Data
	names(lis) <- names(from)
	lis
})
setMethod("[", signature(x = "inverseCDFs", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@.Data <- as(x, "list")[i]
  x@size <- if(is.nothing(x@size)) x@size else x@size[i]
  x@mass0 <- if(is.nothing(x@mass0)) x@mass0 else x@mass0[i]
  stopifnot(validObject(x))
  x
})
setMethod("plot", signature(x = "inverseCDFs", y = "missing"), function(x, y, file = paste("inverseCDFs", Sys.Date(), max(x@size), "pdf", sep = "."), call.par = FALSE, ...)
{
	pdf(file = file)
	if(call.par)
		par(mfrow = c(2, 2))
	for(i in 1:length(x))
	{
		icdf <- x[[i]]
		if(!is.nothing(icdf))
		{
			size <- x@size[i]
			main <- paste(paste(icdf@type, collapse = "; "), size)
			message("\n", i, " ", main, "; file: ", file)
			plot(icdf, main = main, sub = file, ...)
		}
	}
	dev.off()
})
setMethod("Plot", signature(x = "inverseCDFs", y = "missing"), function(x, y, call.par = FALSE, xlab = "estimated local FDR", ylab = "posterior median (log2)", add = FALSE, ...)
{
	if(call.par)
		par(mfrow = c(2, 2))
	med <- Median(x, base = 2)
	graph <- function(fun, ...){fun(x@mass0, med, ...)}
	if(add)
		graph(fun = points, ...)
	else
		graph(fun = plot, xlab = xlab, ylab = ylab, ...)
})
setMethod("Plot", signature(x = "inverseCDFs", y = "xprnSetObject"), function(x, y, xlab = "estimated local FDR", ylab = "estimated mean (log2)", ...)
{
	me <- new.log(stat(y), base = 2)
	nam <- intersect(names(x), names(y))
	plot(x@mass0[nam], me[nam], xlab = xlab, ylab = ylab, col = "gray", ...)
	Plot(x, add = TRUE, ...)
})

setMethod("median", signature(x = "inverseCDFs", na.rm = "missing"), function(x, na.rm)
{
	med <- sapply(x, function(fun){fun(0.5)})
	names(med) <- names(x)
	med
})
setMethod("Median", signature(x = "inverseCDFs"), function(x, base = default(2, "base"), ...)
{
	med <- new.log(median(x = x, ...), base = base)
	med
})



# near end of file:

#Source(file = "interval.s") # functions











