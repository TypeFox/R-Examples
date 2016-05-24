# estimate.r created by David Bickel on 30 August 2007.
# Statomica beats STATOME, STATOMIX, StatomicR, STATOMICS.r, StatomicsPLUS, reSTATOMICS.
# setMethod("annotation", signature(object = "unsortableEstimate") changed on 1 April 2008, after file given to Zahra Montazeri.

#Source("data.r") # for oneLeftOut
setClass("Estimator", representation("function", functional = "logical", type = "character", annotation = "character", nsamples = "Numeric"))
setValidity("Estimator", function(object)
{
	sam.ok <- length(object@nsamples) >= 1
	len.ok <- length(object@functional) == 1 && length(object@type) == 1 && length(object@annotation) == 1
	type.ok <- len.ok && (object@type %in% EstimatorTypes)
	oks <- c(len.ok, type.ok, sam.ok)
	ok <- all(oks==T)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setClass("unsortableEstimate", representation("numeric", annotation = "character", estimator = "Estimator"))
setValidity("unsortableEstimate", function(object)
{
	len.ok <- length(object@annotation) <= 1
	nam <- names(object)
	nam.ok <- is.character(nam)
	if(!nam.ok) warning(class(object), " object has no names, possibly due to creation of sortableEstimate")
	ok <- len.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setClass("estimate", representation(Location = "unsortableEstimate", Scale = "unsortableEstimate", annotation = "character")) #, Size = "Numeric"))
setValidity("estimate", function(object)
{
	vec.names <- setdiff(slotNames(object), "annotation")
	len.ok <- all(length(object@Location) == sapply(vec.names, function(name){length(slot(object, name))})) && length(object@annotation) <= 1
	nam <- names(object)
	nam.ok <- len.ok && is.character(nam) && all(sapply(vec.names, function(name){all(nam == names(slot(object, name)))}))
	ok <- nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("names", signature(x = "estimate"), function(x)
{
	names(Location(x))
})
setMethod("[", signature(x = "estimate", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@Location <- x@Location[i]
	x@Scale<- x@Scale[i]
	x
})
setMethod("sort", signature(x = "estimate"), function(x, decreasing = logical(0), ...)
{#XXX|:removed decreasing = "ANY"
	if(missing(decreasing) || length(decreasing) == 0)
		decreasing <- FALSE
	x[Order(object = x, decreasing = decreasing, ...)]
})

EstimatorTypes <- c("scale", "relative.scale", "inverse.relative.scale", "variance", "relative.variance", "location", "association", "distribution") # modified 24 April 2008 and 20 November 2008

setMethod("annotation", signature(object = "Estimator"), function(object){object@annotation})
zeroHat <- new("Estimator", function(x, y, ...)
{
	0
}, functional = TRUE, type = "location", annotation = "0", nsamples = Numeric(c(1, 2)))
meanHat <- new("Estimator", function(x, y, ...)
{
	sample.mean <- function(object) {mean(object, na.rm = TRUE, ...)}
	if(missing(y))
		sample.mean(x)
	else
		sample.mean(x) - sample.mean(y)
}, functional = TRUE, type = "location", annotation = "sample mean", nsamples = Numeric(c(1, 2)))
meanHat0 <- new("Estimator", function(x, y, ...)
{
	0
}, functional = TRUE, type = "location", annotation = "zero mean", nsamples = Numeric(c(1, 2)))
relative.frequencyHat <- new("Estimator", function(...){relative.frequency(..., threshold = 0, offset.p = NULL)}, functional = TRUE, type = "location", annotation = "relative frequency", nsamples = Numeric(1))
relative.frequencyOffsetHat <- new("Estimator", function(...){relative.frequency(..., threshold = 0, offset.p = 0.5)}, functional = TRUE, type = "location", annotation = "relative frequency (adjusted)", nsamples = Numeric(1)) # relative.frequencyOffsetHat added 090401
aucHat <- new("Estimator", function(...){auc(...)}, functional = TRUE, type = "location", annotation = "P(x > y)", nsamples = Numeric(2))
SdHat <- new("Estimator", function(x, y, var.equal, ...) # estimates sd of x + y or x - y
{
	get.sd <- function(object){Sd(object, ...)}
	if(missing(y))
		get.sd(x)
	else
	{
		if(missing(var.equal))
			var.equal <- default(FALSE, "SdHat's var.equal")
		if(var.equal)
		{
			xy <- c(x - mean(x), y - mean(y))
			sqrt(2) * get.sd(xy)
		}
		else
			sqrt(get.sd(x) ^ 2 + get.sd(y) ^ 2)
	}
}, functional = FALSE, type = "scale", annotation = "st. dev. of diff.", nsamples = Numeric(c(1, 2)))
SeMeanHat <- new("Estimator", function(...){Sd(..., se.of.mean = TRUE)}, functional = FALSE, type = "scale", annotation = "SE of mean", nsamples = Numeric(1))
cvHat <- new("Estimator", function(...){cv(...)}, functional = FALSE, type = "relative.scale", annotation = "coefficient of variation", nsamples = Numeric(1))
functional.cvHat <- new("Estimator", function(...){functional.cv(...)}, functional = TRUE, type = "relative.scale", annotation = "functional coefficient of variation", nsamples = Numeric(1))
t.statHat <- new("Estimator", function(...){t.stat(...)}, functional = FALSE, type = "inverse.relative.scale", annotation = "t statistic", nsamples = Numeric(1))
varHat <- new("Estimator", function(...){var(..., na.rm = TRUE)}, functional = FALSE, type = "variance", annotation = "unbiased variance", nsamples = Numeric(1))
relativeVarHat <- new("Estimator", function(...){relativeVar(...)}, functional = FALSE, type = "relative.variance", annotation = "relative variance", nsamples = Numeric(1))
functional.relativeVarHat <- new("Estimator", function(...){functional.relativeVar(...)}, functional = TRUE, type = "relative.variance", annotation = "functional relative variance", nsamples = Numeric(1))
relative.frequencyHat.norm <- new("Estimator", function(x, mu, stdev, ...) # cf. relative.frequencyHat
{
	if(missing(mu)) mu <- 0
	sample.mean <- function(object) {mean(object, na.rm = TRUE, ...)}
	mean <- sample.mean(x)
	if(missing(stdev))
		stdev <- try(sd(x, na.rm = TRUE))
	if(is(stdev, "try-error"))
	{ message("bad stdev in relative.frequencyHat.norm"); browser()}
	pnorm(q = mu, mean = mean, sd = stdev, lower.tail = FALSE)
}, functional = TRUE, type = "location", annotation = "fraction greater", nsamples = Numeric(1))
aucHat.norm0 <- aucHat0 <- relative.frequencyHat.norm0 <- relative.frequencyHat0 <- new("Estimator", function(x, y, ...)
{
	0.5
}, functional = TRUE, type = "location", annotation = "fraction greater is 50%", nsamples = Numeric(c(1, 2)))
aucHat.norm <- new("Estimator", function(x, y, mu, stdev, var.equal, ...) # tested in fma090115.r
{
	if(missing(mu)) mu <- 0
	sample.mean <- function(object) {mean(object, na.rm = TRUE, ...)}
	mean <- sample.mean(x) - sample.mean(y)
	if(missing(stdev))
	{
		if(missing(var.equal))
			var.equal <- default(FALSE, "aucHat.norm's var.equal")
		stdev <- try(SdHat(x = x, y = y, var.equal = var.equal, na.rm = TRUE))
	}
	if(is(stdev, "try-error"))
	{ message("bad stdev in aucHat.norm"); browser()}
	pnorm(q = mu, mean = mean, sd = stdev, lower.tail = FALSE)
}, functional = TRUE, type = "location", annotation = "fraction greater", nsamples = Numeric(2))
aucHat.norm.vareq <- new("Estimator", function(...)
{
	aucHat.norm(..., var.equal = TRUE)
}, functional = TRUE, type = "location", annotation = "fraction greater", nsamples = Numeric(2))


#setClass("unsortableEstimate", representation("numeric", annotation = "character", estimator = "Estimator"))
#setValidity("unsortableEstimate", function(object)
#{
#	len.ok <- length(object@annotation) <= 1
#	nam <- names(object)
#	nam.ok <- is.character(nam)
#	if(!nam.ok) warning(class(object), " object has no names, possibly due to creation of sortableEstimate")
#	ok <- len.ok
#	if(!ok)
#	{ printInvalid(object); browser()}
#	ok
#})
setMethod("annotation", signature(object = "unsortableEstimate"), function(object)
{
	if(length(object@annotation) == 0)
		object@estimator@annotation
	else
		paste(object@annotation, object@estimator@annotation, sep = "; ")
})
setAs(from = "unsortableEstimate", to = "numeric", function(from)
{
	vec <- as.numeric(from)
	names(vec) <- names(from)
	vec
})
setMethod("[", signature(x = "unsortableEstimate", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	vec <- as(x, "numeric")[i]
#	names(vec) <- names(x)[i]
	if(is.null(names(vec)))
	{ message("'[' names error"); browser()}
	new.unsortableEstimate(vec, annotation = x@annotation, estimator = x@estimator)
})

setClass("sortableEstimate", representation("unsortableEstimate", LocationScale = "estimate"))
setValidity("sortableEstimate", function(object)
{
	len.ok <- length(object) == length(object@LocationScale@Scale)
	nam.ok <- len.ok && all(names(object) == names(object@LocationScale))
	ok <- nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setAs(from = "sortableEstimate", to = "unsortableEstimate", function(from)
{
	vec <- as.numeric(from)
	if(is.null(names(vec)))
	{
		names(vec) <- if(is.null(names(from))) # this occurs with new("unsortableEstimate", ...)
			names(from@LocationScale)
		else
			names(from)
	}
	if(is.null(names(vec)))
	{ warning("from = 'sortableEstimate' vec lacks names")}
	new.unsortableEstimate(vec, annotation = from@annotation, estimator = from@estimator)
})
setAs(from = "sortableEstimate", to = "numeric", function(from)
{
	as(as(from, "unsortableEstimate"), "numeric")
})
Extract.sortableEstimate <- function(x, i, j, drop)
{
	stopifnot(is(x, "sortableEstimate"))
	if(is.null(names(x)))
	{ message("'[' sortableEstimate names error"); browser()}
	vec <- as(x, "unsortableEstimate")[i]
	LocationScale <- x@LocationScale[i]
	new.sortableEstimate(vec, LocationScale = LocationScale)
}
setMethod("[", signature(x = "sortableEstimate", i = "ANY", j = "missing"), Extract.sortableEstimate)
setMethod("sort", signature(x = "sortableEstimate"), function(x, decreasing = logical(0), ...)
{#XXX|:removed decreasing = "ANY"
	if(missing(decreasing) || length(decreasing) == 0)
		decreasing <- FALSE
	if(is.null(names(x)))
	{ message("sort sortableEstimate names error"); browser()}
	x[Order(object = x, decreasing = decreasing, ...)]
})

setClassUnion("numericEstimate", c("numeric", "unsortableEstimate", "sortableEstimate")) # used in slot of biasEstimate class

setClass("EstimateLeftOut", representation("oneLeftOut", shrinkage = "Scalar"))
setValidity("EstimateLeftOut", function(object)
{
	is(noneLeftOut(object), "unsortableEstimate")
})

setClass("biasEstimates", representation("list", shrinkage = "numeric"))
setValidity("biasEstimates", function(object)
{
	cla.ok <- all(sapply(object, is, class2 = "biasEstimate"))
	len.ok <- length(object@shrinkage) == length(object)
	a.ok <- all(annotation(object) == sapply(object, annotation))
	ok <- cla.ok && len.ok && a.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("annotation", signature(object = "biasEstimates"), function(object){annotation(object[[1]])})
setMethod("plot", signature(x = "biasEstimates", y = "missing"), function(x, y, col, ranks, nfeatures, use.lower.ranks, file, call.biasEstimates, smooth, main, smooth.globally, f, save.space, cex, call.par, delay.legend, ...)
{
	if(missing(call.par)) call.par <- TRUE
	if(missing(save.space))
	{
		save.space <- FALSE
		message("save.space = ", save.space)
	}
	if(missing(delay.legend))
	{
		delay.legend <- FALSE
		if(save.space)
			message("delay.legend = ", delay.legend)
	}
	if(missing(cex))
	{
		cex <- if(save.space)
			1.5
		else
			0.95 # par()$cex * 1
	}
	smooth.x <- function(smooth, x)
	{
		if(is.logical(smooth) && !smooth)
			x # stop("cannot smooth x")
		else
		{
			if(!is.function(smooth))
				smooth <- smoothDataFUN(f = f)
			smooth(x)
		}
	}
	if(missing(smooth.globally))
	{
		smooth.globally <- missing(smooth) || !is.logical(smooth) || smooth
		message("smooth.globally = ", smooth.globally)
	}
	if(missing(f))
	{
		f <- if(smooth.globally) 1/100 else 1
		if(missing(smooth) || !identical(smooth, FALSE))
			message("f = ", f)
	}
	if(missing(main))
		main <- annotation(x)
	if(missing(call.biasEstimates))
	{
		call.biasEstimates <- TRUE
		message("call.biasEstimates = ", call.biasEstimates)
	}
	pch <- 1:length(x@shrinkage)
	if(missing(col))
		col <- pch
	if(length(col) == 1)
		col <- rep(col, length(x@shrinkage))
	stopifnot(length(col) == length(pch))
	if(call.par)
		par(mfrow = c(2, 2))
	total.nfeatures <- length(x[[1]])
	if(missing(ranks))
	{
		if(missing(nfeatures) && call.biasEstimates)
		{
			nfeatures <- min(200, total.nfeatures)
			message("nfeatures = ", nfeatures)
		}
		graph <- function(use.lower.ranks, smooth, smooth.globally, call.par)
		{
			if(smooth.globally)
			{
				graph.object <- smooth.x(smooth = smooth, x = x)
				smooth <- FALSE
			}
			else
				graph.object <- x
			plot(x = biasEstimates(object = graph.object, nfeatures = ceiling(nfeatures / 2), use.lower.ranks = use.lower.ranks), col = col, file = character(0), call.biasEstimates = FALSE, smooth = smooth, smooth.globally = smooth.globally, main = main, f = f, save.space = save.space, cex = cex, call.par = call.par, delay.legend = delay.legend, ...)
		}
		if(missing(file) && call.biasEstimates && save.space)
		{
			if(!missing(use.lower.ranks))
				message("use.lower.ranks ignored")
			if(missing(smooth)) smooth <- TRUE
			graph(TRUE, smooth = smooth, smooth.globally = smooth.globally, call.par = call.par)
			graph(FALSE, smooth = smooth, smooth.globally = smooth.globally, call.par = FALSE)
			return(character(0))
		}
		else if(is.character(file) && length(file) == 1)
		{
			if(!call.biasEstimates)
				stop("file specified, but call.biasEstimates is not TRUE")
			if(!missing(use.lower.ranks))
				message("use.lower.ranks ignored")
			pdf(file = file)
			if(missing(smooth)) smooth <- TRUE
			graph(TRUE, smooth = FALSE, smooth.globally = FALSE, call.par = call.par)
			graph(FALSE, smooth = FALSE, smooth.globally = FALSE, call.par = FALSE)
			graph(TRUE, smooth = smooth, smooth.globally = smooth.globally, call.par = FALSE)
			graph(FALSE, smooth = smooth, smooth.globally = smooth.globally, call.par = FALSE)
			dev.off()
			return(file)
		}
		else
			message("input parameters will not yield the ideal plot")
#		else if(is.ch
#		{
#			message("bad file argument:")
#			print(file)
#			browser()
#		}
		if(missing(use.lower.ranks) && call.biasEstimates)
		{
			use.lower.ranks <- logical(0)
			message("plot's use.lower.ranks = ", use.lower.ranks)
		}
		if(smooth.globally)
		{
			x <- smooth.x(smooth = if(missing(smooth)) TRUE else smooth, x = x)
			smooth <- FALSE
		}
		if(call.biasEstimates)
			x <- biasEstimates(x, nfeatures = nfeatures, use.lower.ranks = use.lower.ranks)
		ranks <- Rank(x)
	}
	else # ranks specified
	{
		if(!missing(file))
			warning("file not used since ranks was specified")
		if(smooth.globally)
		{
			x <- smooth.x(smooth = if(missing(smooth)) TRUE else smooth, x = x)
			smooth <- FALSE
		}
		if(!is.numeric(ranks))
		{
			ranks <- Rank(x)
			message("ranks:")
			print(summary(ranks))
		}
		if(!all(ranks == Rank(x)) && call.biasEstimates)
			x <- biasEstimates(x, ranks = ranks)
	}
	if(!is(x, "biasEstimates")) stop("x lost its class")
	if(!all(ranks == Rank(x)))
		stop("plot biasEstimates error")
	type <- if(!missing(smooth) && (!is.logical(smooth) || smooth))
	{
		x <- smooth.x(smooth = smooth, x = x)
		"l"
	}
	else
	  "p"
	if(smooth.globally)
		type <- "l"
	if(!is(x, "biasEstimates")) stop("smooth x lost its class")
	if(!(type %in% c("p", "l"))) stop("bad type")
	supergraph <- function(x.fun, y.fun, xlab, ylab, type)
	{
		get.vec <- function(be, fun)
		{
			stopifnot(is(be, "biasEstimate"))
			vec <- as.numeric(fun(be))
#			vec <- biasEstimate(object = x, ranks = ranks) # vec[ranks]
			stopifnot(length(vec) == length(ranks))
			vec
		}
		subgraph <- function(gfun, x, y, col, pch)
		{
			sub <- function(..., pch)
			{
				if(type == "p")
					gfun(..., pch = pch)
				else if(type == "l")
					gfun(..., lty = pch)
				else
					stop("very bad type 2")
			}
			sub(x = x, y = y, col = col, pch = pch, ...)
		}
		if(missing(x.fun) && missing(xlab))
			xlab <- "feature rank"
		else if(is.function(x.fun) && is.character(xlab))
		{
			xs <- try(sapply(x, get.vec, fun = x.fun))
			if(is(xs, "try-error"))
			{ message("bad xs"); browser()}
		}
		else
			stop("bad x.fun and/or xlab")
		ys <- sapply(x, get.vec, fun = y.fun)
		get.lim <- function(mat)
		{
			vec <- as.numeric(mat)
			range(vec[is.finite(vec)], na.rm = TRUE)
		}
		if(missing(x.fun))
		{
			xlim <- range(ranks)
			ylim <- get.lim(ys)
		}
		else
		{
			xlim <- ylim <- get.lim(rbind(xs, ys))
		}
		for(j in 1:length(x))
		{
			gfun <- if(j == 1) # || !missing(x.fun))
				function(...){plot(..., xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = ylim, type = type)}
			else if(type == "p")
				points
			else if(type == "l")
				lines
			else
				stop("very bad type 1")
			xvec <- if(missing(x.fun))
				ranks
			else
			  xs[, j]
			subgraph(gfun = gfun, x = xvec, y = ys[, j], col = col[j], pch = pch[j])
			if(!missing(x.fun))
				abline(a = 0, b = 1, col = "yellow")
		}
		as.numeric(ys) # used to set legend.x
	}
	add.legend <- function(legend.x)
	{
		leg <- function(...)
		{
			if(save.space)
				legend(..., cex = cex)
			else if(type == "p")
				legend(..., pch = pch)
			else if(type == "l")
				legend(..., lty = pch, cex = cex)
			else
				stop("very bad type 3")
		}
		rounded.shrinkage <- round(x@shrinkage, 3)
		if(save.space)
		{
			get.expression <- function(r.s)
			{
				cmd <- paste("expression(paste(italic(zeta), '=', ", r.s, "))")
				eval(parse(text = cmd))
			}
			leg(legend = sapply(rounded.shrinkage, get.expression), x = legend.x, text.col = col, ncol = 3, bty = "n") # title = expression(zeta) horiz = TRUE
		}
		else
			leg(legend = paste("shrinkage parameter =", rounded.shrinkage), x = legend.x, col = col, bty = "n")
	}
	if(!save.space)
		supergraph(x.fun = function(object){object@uncorrected}, y.fun = corrected, xlab = "uncorrected estimate", ylab = "corrected estimate", type = "p")
	ys <- supergraph(y.fun = corrected, ylab = "corrected estimate", type = "p")
	add.saving.legend <- function()
	{
		midpoint <- mean(range(ys, na.rm = TRUE))
		midpoint.excess <- midpoint - median(ys, na.rm = TRUE)
		message("midpoint exceeds median by ", midpoint.excess)
		add.legend(legend.x = if(midpoint.excess > 0) "top" else "bottom")
	}
	if(save.space && !delay.legend)
		add.saving.legend()
	supergraph(y.fun = identity, ylab = "estimated bias", type = type)
	if(save.space && delay.legend)
		add.saving.legend()
	if(!save.space)
	{
		blank.plot()
		add.legend(legend.x = "center")
	}
})

setClass("sim.biasEstimates", representation(gene.rank = "Numeric", uncorrected.estimate = "numeric", rough.bias.estimate = "numeric", smooth.bias.estimate = "numeric", sample.size = "Scalar", means = "numeric", sds = "Numeric", rxprnSet = "function", estimator = "Estimator", shrinkage = "Scalar", true.rank = "Scalar")) # for a single gene and many samples simulated according to rxprnSet
setValidity("sim.biasEstimates", function(object)
{
	len.ok <- all(length(object@gene.rank) == sapply(c("uncorrected.estimate", "rough.bias.estimate", "smooth.bias.estimate"), function(name){length(slot(object, name = name))}))
	ngenes.ok <- TRUE # length(object@means) == length(object@sds)
	fun.ok <- is(object@rxprnSet(), "xprnSet")
	ok <- len.ok && ngenes.ok && fun.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("print", signature(x = "sim.biasEstimates"), function(x)
{
	stop("implement me")
})
setMethod("[", signature(x = "sim.biasEstimates", i = "ANY", j = "ANY", drop = "missing"), function(x, i, j, drop)
{
	new.sim.biasEstimates(gene.rank = x@gene.rank[i], uncorrected.estimate = x@uncorrected.estimate[i], rough.bias.estimate = x@rough.bias.estimate[i], smooth.bias.estimate = x@smooth.bias.estimate[i], sample.size = x@sample.size, means = x@means[j], sds = x@sds[j], rxprnSet = x@rxprnSet, estimator = x@estimator, shrinkage = x@shrinkage, true.rank = x@true.rank)
})
setMethod("hist", signature(x = "sim.biasEstimates"), function(x, call.par, breaks, length.out, include.mean, ...)
{
	if(missing(include.mean))
	{
		include.mean <- FALSE
		message("include.mean = ", include.mean)
	}
	if(missing(call.par))
		call.par <- TRUE
	if(call.par)
		par(mfrow = c(2, 2))
	if(missing(breaks))
		breaks <- NULL
	if(missing(length.out))
	{
		length.out <- 20
		message("length.out = ", length.out)
	}
	get.breaks <- function(x)
	{
		num.breaks <- length.out
		seq(min(x), max(x), length.out = num.breaks)
	}
	get.hist <- function(x, ...)
	{
		if(is.null(breaks))
			breaks <- get.breaks(x = x)
		hist(x = x, main = "histogram", ylab = "number of genes", breaks = breaks, ...)
	}
	tm.lab <- "true mean of log(expression ratio)"
	tdel.lab <- "true differential expression level"
	if(include.mean)
		get.hist(x = x@means, xlab = tm.lab, ...)
	tdel <- estimand.value(x, gene = 1:length(x@means))
	get.hist(x = tdel, xlab = tdel.lab, ...)
	get.ecdf <- function(x, ...)
	{
		obj <- ecdf(x)
		plot(obj, ylab = "proportion of genes with lower levels", main = "empirical distribution function", ...)
	}
	if(include.mean)
		get.ecdf(x = x@means, xlab = tm.lab)
	get.ecdf(x = tdel, xlab = tdel.lab)
})
setMethod("plot", signature(x = "sim.biasEstimates", y = "missing"), function(x, y, bias, biased.estimate, call.par, xlim, ...)
{
	if(missing(call.par))
		call.par <- TRUE
	if(call.par)
		par(mfrow = c(2, 2))
	true.level <- estimand.value(object = x)
	if(missing(bias) && missing(biased.estimate))
	{
		bias <- true.bias(x)
		message("bias = ", bias)
	}
	else if(missing(bias))
	{
		bias <- biased.estimate - true.level
	}
	else
		stopifnot(missing(biased.estimate))
	if(missing(biased.estimate))
		biased.estimate <- true.level + bias
	get.main <- function(digits)
	{
		ro <- function(object){round(object, digits = digits)}
		paste("true bias = ", ro(biased.estimate), " n ", ro(true.level), " = ", ro(bias), sep = "")
	}
	conditional.bias.plot <- function(xlim)
	{
		plot(x = x@gene.rank, y = x@uncorrected.estimate, xlab = "estimated rank", ylab = "uncorrected estimate", xlim = xlim, main = get.main(2))
		abline(v = x@true.rank)
		abline(h = biased.estimate, col = "gray")
		abline(h = true.level, col = "black")
	}
	conditional.bias.plot(xlim = range(c(x@gene.rank, x@true.rank)))
	if(missing(xlim))
		message("specify xlim for another conditional bias plot")
	else
		conditional.bias.plot(xlim = xlim)
	get.lab <- function(mse, prefix)
	{
		paste(prefix, " bias estimate (RMSE: ", round(sqrt(mse), 3), ")", sep = "")
	}
	get.mse <- function(vec)
	{
		mean.squaredError(predicted = vec, observed = rep(bias, length(vec)))
	}
	rough <- x@rough.bias.estimate
	smoothe <- x@smooth.bias.estimate
	rough.mse <- get.mse(rough)
	smooth.mse <- get.mse(smoothe)
	xlim <- ylim <- range(c(rough, smoothe), na.rm = TRUE)
	main <- paste("MSE reduction: ", 100 * round((rough.mse - smooth.mse) / rough.mse, 3), "%", sep = "")
	plot(x = rough, y = smoothe, xlab = get.lab(rough.mse, "rough"), ylab = get.lab(smooth.mse, "smooth"), xlim = xlim, ylim = ylim, main = main, ...)
	bias.col <- "gray"
	abline(v = bias, col = bias.col)
	abline(h = bias, col = bias.col)
#	abline(a = 0, b = 1, col = bias.col)
})

# used by extended.r and/or extended.s:
setMethod("hist", signature(x = "xprnSet"), function(x, ...)
{
	hist(x = mean(x), main = annotation(x), ...)
})

setMethod("plot", signature(x = "xprnSet", y = "missing"), function(x, ...)
{
	plot(x = Density(mean(x)), ...)
})


# near end of file:

#Source(file = "estimate.s") # functions
