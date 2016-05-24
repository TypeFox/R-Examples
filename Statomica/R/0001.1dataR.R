# Created by David Bickel on 17 August 2007.
# setMethod("[", signature(x = "Matrix"... changed on 1 Feb. 2008.
# setReplaceMethod("[", "Matrix"... added on 3 April 2008.
# setReplaceMethod("[", "Numeric"... added on 3 April 2008.
# "plot", signature(x = "numeric", y = "biasEstimate") changed 9 April 2008.
# More changes 16 April 2008 and thereafter.
# April 2013: marta added featureData to inputs of xprnSet and XprnSet to be used with reference sets
#library(Biobase)
#library(multtest)
#library(graphics)
#library(base)
removeMethods.boo<-FALSE
if(!is.logical(try(removeMethods.boo)))
	removeMethods.boo <- FALSE#TRUE
if(removeMethods.boo)
{
	message("removing methods; set removeMethods.boo <- FALSE to disable")
	try(removeMethods("plot")) # added 2 July 2008
	try(removeMethods("stripchart")) # added 20 Aug. 2008
	try(removeMethods("summary")) # added 8 January 2009
	try(removeMethods("print"))
	try(removeMethods("mean"))
	try(removeMethods("exp"))
}



hvline <- function(x, y, h, save.time = default(FALSE, "save.time"), ...)
{
	col <- "gray"
	if(save.time)
		abline(h = h, col = col)
	else if(length(h) > 0 && is.numeric(h) && length(y) == length(x))
	{
		for(H in h)
		{
			diffs <- abs(y - H)
			v <- x[diffs == min(diffs)]
			v <- v[v != min(x) & v != max(x)]
			if(length(v) > 0)
			{
				cat(class(x), ": drawing vertical line at "); print(v)
				abline(h = H, v = v, col = col)
			}
		}
	}
	else
		message("cannot draw lines")
}

setMethod("plot", signature(x = "numeric", y = "function"), function(x, y, type, add = default(FALSE, "add"), call.browser = FALSE, h, save.time = default(FALSE, "save.time"), ylim, ...)
{
	x <- sort(x)
	if("ylim" %in% names(list(...))) { message("ylim in ... (", class(x), ", ", class(y), ")"); browser()}
	if(missing(h))
		h <- default(if(is(y, "CDF")) c(0.25, 0.5, 0.75) else numeric(0), paste(class(y), '"plot", signature(x = "numeric", y = "function")'))
	if(missing(type))
		type <- default("l", "type")
	if(call.browser)
	{ message("press c to continue; type: ", type); browser()}
	Y <- sapply(x, y)
	if(missing(ylim) || is.null(ylim))
		ylim <- range(Y, na.rm = TRUE)
	fun <- if(add)
	{
		if(type == "l")
			lines
		else if(type == "p")
			points
		else
			stop('bad args in "plot", signature(x = "numeric", y = "function")')
	}
	else
		function(...){plot(..., ylim = ylim)}
	message("plotting Y of class ", class(x), "\n")
  fun(x = x, y = Y, type = type, ...)
  hvline(x = x, y = Y, h = h, save.time = save.time)
  invisible(Y)
})
setMethod("plot", signature(x = "numeric", y = "list"), function(x, y, ...) # usage: plot(x = c(1,2), y = list(a = c(3,4), b = c(5,6)), legend.x = "topleft")
{
	get.col <- function(elem)
	{
		vec <- if(is.function(elem))
			elem(x)
		else
			as(elem, "numeric")
		stopifnot(length(x) == length(vec))
		vec
	}
	plot(x = x, y = data.frame(lapply(y, get.col)), ...) # lapply(y, as, Class = "numeric"))
})
setMethod("plot", signature(x = "missing", y = "numeric"), function(x, y, sub, ...)
{
	qcol <- "gray"
	if(missing(sub))
		sub <- paste("(", qcol, " lines are quantiles)", sep = "")
	x <- 1:length(y)
	plot(x = x, y = y, sub = sub, xlab = "index", xlim = range(x), ...)
	abline(h = quantile(y), col = qcol)
	abline(h = median(y), col = gray(level = 0.4))
})
setMethod("plot", signature(x = "numeric", y = "data.frame"), function(x, y, col, pch, type, legend.legend, legend.x, call.legend = default(ncol(y) > 1, "call.legend"), xlim, ylim, ylab = default(if(ncol(y) == 1) names(y) else "y"), normalizeFUN, transformFUN = NULL, MARGIN, log = "", text.offset, countable.values = NULL, cex = default(0.7, "cex"), multiple = FALSE, abline.args = list(), call.par = multiple, main, ...)
{
	if("ylim" %in% names(list(...))) { message("ylim in ... (", class(x), ", ", class(y), ")"); browser()}
	stopifnot(length(x) == nrow(y))
	if(missing(MARGIN))
		MARGIN <- 1
	if(is.function(transformFUN))
	{
		for(j in 1:ncol(y))
			y[, j] <- transformFUN(y[, j])
		message("data frame transformed")
	}
	if(!missing(normalizeFUN) && is.function(normalizeFUN))
	{
		if(MARGIN == 1)
		{
			for(i in 1:nrow(y))
				y[i, ] <- y[i, ] / normalizeFUN(y[i, ])
		}
		else if(MARGIN == 2)
		{
			for(j in 1:ncol(y))
				y[, j] <- y[, j] / normalizeFUN(y[, j])
		}
		else
			stop("bad MARGIN")
	}
	if(missing(col))
		col <- 1 # :length(x)
	if(missing(pch))
		pch <- 1:length(x)
	if(missing(type))
		type <- default("p", "type")
	if(missing(legend.x))
		legend.x <- default("center", "legend.x")
	if(length(col) == 1)
		col <- rep(col, length(x))
	if(length(pch) == 1)
		pch <- rep(pch, length(x))
	nam <- colnames(y)
	if(missing(legend.legend))
		legend.legend <- default(if(is.null(nam)) col else nam, "legend.legend")
	get.lim <- function(vec, is.positive)
	{
		vec <- as(vec, "numeric")
		assert.is(vec, "numeric")
		get.vec <- function(vec, boo)
		{
			boo <- ifelse(is.na(boo), FALSE, boo)
			if(all(boo))
				vec
			else
			{
				message("\nAt least ", sum(!boo), " values were removed from compuation of xlim or ylim:")
				pri <- try(print(stats(vec[!boo])))
				if(is.err(pri))
				{ message("cannot print removed values")}
				cat("\n")
				ifelse(boo, vec, min(vec[boo]))
			}
		}
		boo <- if(is.positive)
			is.finite(logb(vec))
		else
			is.finite(vec)
		if(is.positive && identical(transformFUN, exp) && !all(boo, na.rm = TRUE))
		{ message("transformFUN error"); browser()}
		vec <- get.vec(vec = vec, boo = boo)
		lim <- range(vec, na.rm = TRUE)
		stopifnot(length(lim) == 2)
		if(is.positive)
			stopifnot(all(lim > 0))
		lim
	}
	if(missing(xlim))
		xlim <- get.lim(x, is.positive = log %in% c("x", "xy"))
	stopifnot(length(x) == nrow(y))
	if(missing(ylim))
		ylim <- get.lim(as.numeric(as.matrix(y[x >= xlim[1] & x <= xlim[2], ])), is.positive = log %in% c("y", "xy"))
	if(call.par)
		Mfrow()
	graph <- function(...)
	{
		plot(..., type = type, xlim = xlim, ylim = ylim, ylab = ylab, log = log)
		if(length(abline.args) > 0)
			do.call(abline, abline.args)
	}
	plo <- try(graph(x = x, y = y[[1]], col = col[1], pch = pch[1], ...))
	if(is.err(plo))
	{ message("bad plo"); browser()}
	fun <- if(multiple)
	{
		title(main = if(missing(main)) legend.legend[1] else main)
		col <- rep(col[1], length(col))
		pch <- rep(pch[1], length(pch))
		graph
	}
	else if(type == "p")
		points
	else if(type == "l")
		lines
	else
		stop("bad type")
	if(ncol(y) > 1)
	{
		for(j in 2:ncol(y))
		{
			yvec <- y[, j]
			fun(x = x, y = yvec, col = col[j], pch = pch[j], ...)
			title(main = if(missing(main)) legend.legend[j] else main)
		}
	}
	if(!missing(text.offset))
	{
		stopifnot(length(text.offset) == 1)
		write.text <- function(vec, ax)
		{
			un <- unique(vec)
			if(!is.null(countable.values))
			{
				stopifnot(is.numeric(countable.values) && length(countable.values) == 1)
				un <- un[un %in% countable.values]
			}
			if(length(un) >= 1)
			{
				offsets <- rep(text.offset, length(un))
				text.fun <- function(X, Y)
				{
					X <- as(X, "numeric")
					Y <- as(Y, "numeric")
					labels <- sapply(un, function(un1){sum(un1 == vec)})
					datf <- try(data.frame(x = X, y = Y, labels = labels))
					if(!is(datf, "data.frame"))
					{ message("bad datf"); browser()}
					print(datf)
					text(x = X, y = Y, labels = labels, cex = cex)
				}
				if(ax == "x")
					text.fun(X = un, Y = offsets)
				else if(ax == "y")
					text.fun(X = offsets, Y = un)
				else
					stop("bad ax")
			}
		}
		write.text(vec = x, ax = "x")
		j <- 1
		yvec <- y[, j]
		write.text(vec = yvec, ax = "y")
		if(ncol(y) > 1)
			message("only numbers of element ", j, " of y plotted")
	}
	if(call.legend && !multiple && length(legend.x) > 0)
		legend(x = legend.x, legend = legend.legend, col = col, pch = pch)
})

setClass("pNumeric", representation("numeric"))
setValidity("pNumeric", function(object)
{
	num.ok <- all(is.na(object)) || all(object >= 0 & object <= 1, na.rm = TRUE)
	if(!num.ok)
	{ printInvalid(object); browser()}
	num.ok
})
setAs(from = "pNumeric", to = "numeric", function(from)
{
  x <- from@.Data
  names(x) <- names(from)
  x
})
setMethod("[", signature(x = "pNumeric", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  pNumeric(as(x, "numeric")[i])
})
setReplaceMethod("[", signature(x = "pNumeric", i = "ANY", j = "missing"), function(x, i, j, value) # showMethods("[<-")
{
	vec <- as(x, "numeric")
	tr <- try(vec[i] <- value)
	if(is(tr, "try-error") || length(vec) != length(x))
	{ message("bad replacement of ", class(x)); browser()}
	pNumeric(vec)
})

setClass("Numeric", representation("numeric"))
setValidity("Numeric", function(object)
{
	num.ok <- all(is.na(object)) || all(object >= 0, na.rm = TRUE)
	if(!num.ok)
	{ printInvalid(object); browser()}
	num.ok
})
setAs(from = "Numeric", to = "numeric", function(from)
{
  x <- from@.Data
  names(x) <- names(from)
  x
})
setMethod("[", signature(x = "Numeric", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  Numeric(as(x, "numeric")[i])
})
setReplaceMethod("[", signature(x = "Numeric", i = "ANY", j = "missing"), function(x, i, j, value) # showMethods("[<-")
{
	vec <- as(x, "numeric")
	vec[i] <- value
#	tr <- try(vec[i] <- as(value, "numeric"))
#	if(is(tr, "try-error") || length(vec) != length(x))
#	{ message("bad replacement of ", class(x)); browser()}
	Numeric(vec)
})
#removeMethods("sqrt")
setMethod("sqrt", signature(x = "Numeric"), function(x)
{
	vec <- as(x, "numeric")
	nam <- names(x)
  x@.Data <- sqrt(vec)
  names(x) <- nam
  stopifnot(validObject(x))
  x
})

setClass("scalar", representation("numeric"))
setValidity("scalar", function(object)
{
	length(object) %in% c(0, 1) # changed 16 April 2008
})
setClass("Scalar", representation("scalar"))
setValidity("Scalar", function(object)
{
	length(object) == 0 || is.na(object) || object >= 0 # changed 16 April 2008
})
setIs("Scalar", "Numeric")

setClass("Matrix", representation("matrix"))
setValidity("Matrix", function(object)
{
	all(as.numeric(object) >= 0, na.rm = TRUE)
})
setMethod("[", signature(x = "Matrix", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  mat <- as(x, "matrix")
  mat <- mat[i, , drop = FALSE]
  Matrix(mat)
})
setMethod("[", signature(x = "Matrix", i = "missing", j = "ANY"), function(x, i, j, drop)
{
  mat <- as(x, "matrix")
  mat <- mat[, j, drop = FALSE]
  Matrix(mat)
})
setMethod("[", signature(x = "Matrix", i = "ANY", j = "ANY"), function(x, i, j, drop)
{
  mat <- as(x, "matrix")
  mat <- mat[i, j, drop = FALSE]
  Matrix(mat)
})
setReplaceMethod("[", signature(x = "Matrix", i = "ANY", j = "missing"), function(x, i, j, value) # showMethods("[<-")
{
#	cat("replacing with "); print(value)
	mat <- as(x, "matrix")
	tr <- try(mat[i, ] <- value)
	if(is(tr, "try-error"))
	{ message("bad replacement"); browser()}
	Matrix(mat)
})

setClass("xprnSet", representation(es = "ExpressionSet"))
##setIs("xprnSet", "ExpressionSet")
setValidity("xprnSet", function(object)
{
	all(colnames(exprs(object)) == rownames(pData(object)))
})
setAs(from = "xprnSet", to = "ExpressionSet", function(from)
{
  from@es
})
setAs(from = "xprnSet", to = "numeric", function(from)
{
	as(exprs(from), "numeric")
})
setMethod("names", signature(x = "xprnSet"), function(x)
{
	featureNames(x)
})
setMethod("length", signature(x = "xprnSet"), function(x) # like Size
{
	size <- sapply(1:length(names(x)), function(i){sum(is.finite(exprs(x)[i, ]))})
	names(size) <- names(x)
	size
})
setMethod("print", signature(x = "xprnSet"), function(x)
{
  message(class(x), "; high-throughput data:\n")
  print(as(x, "ExpressionSet"))
})
setMethod("exprs", signature(object = "xprnSet"), function(object)
{
  exprs(as(object, "ExpressionSet"))
})
setMethod("featureNames", signature(object = "xprnSet"), function(object)
{
  featureNames(as(object, "ExpressionSet"))
})
setMethod("pData", signature(object = "xprnSet"), function(object)
{
  pData(as(object, "ExpressionSet"))
})
setMethod("[", signature(x = "xprnSet", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  es <- as(x, "ExpressionSet")
  x@es <- es[i, ]
  x
})
setMethod("[", signature(x = "xprnSet", i = "missing", j = "ANY"), function(x, i, j, drop)
{
  es <- as(x, "ExpressionSet")
  x@es <- es[, j]
  x
})
setMethod("[", signature(x = "xprnSet", i = "ANY", j = "ANY"), function(x, i, j, drop)
{
  es <- as(x, "ExpressionSet")
  x@es <- es[i, j]
  x
})
setMethod("logb", signature(x = "ExpressionSet", base = "missing"), function(x, base)
{
  exprs(x) <- logb(exprs(x))
  x
})
setMethod("logb", signature(x = "xprnSet", base = "missing"), function(x, base)
{
  x@es <- logb(as(x, "ExpressionSet"))
  x
})
setMethod("exp", signature(x = "ExpressionSet"), function(x)
{
  exprs(x) <- exp(exprs(x))
  x
})
setMethod("exp", signature(x = "xprnSet"), function(x)
{
  x@es <- exp(as(x, "ExpressionSet"))
  x
})
# setMethod("plot", signature(x = "xprnSet", y = "missing"), function(x, ...) # see estimated.r and Plot
setMethod("plot", signature(x = "xprnSet", y = "xprnSet"), function(x, y, ...)
{
	stopifnot(all(dim(x) == dim(y)))
	plot(x = as(x, "numeric"), y = as(y, "numeric"), xlab = annotation(x), ylab = annotation(y), ...)
})
setMethod("stripchart", signature(x = "xprnSet"), function(x, file, factor.name, ...)
{
	if(missing(file))
	{
		file <- paste(annotation(x), factor.name, "pdf", sep = ".")
		message("file = ", file)
	}
	if(missing(factor.name))
	{
		stop("not yet implemented for missing factor.name")
		fac <- factor()
		levs <- "level"
		message("factor.name missing")
	}
	else
	{
		fac <- pData(x)[, factor.name]
		stopifnot(is.factor(fac))
		levs <- levels(fac)
	}
	call.pdf <- is.character(file) && length(file) == 1
	if(call.pdf) pdf(file = file)
	par(mfrow = c(2, 2))
	for(feature in featureNames(x))
	{
		y <- exprs(x[feature, ])
		boo <- is.finite(y)
		if(any(boo))
		{
			y <- y[boo]
			fac <- fac[boo]
			if(length(y) != length(fac))
			{ message("bad length of y or fac"); browser()}
			stripchart(y ~ fac, main = feature, xlab = feature, ...)
		}
		else
			message("data not available for ", feature)
	}
	if(call.pdf) dev.off()
})
setMethod("dim", signature(x = "xprnSet"), function(x)
{
	dim(exprs(x))
})
setMethod("dimnames", signature(x = "xprnSet"), function(x)
{
	dimnames(exprs(x))
})

setClass("Bernoulli.normal.rxprnSet", representation("xprnSet", alternative = "logical", param = "numeric"))
setValidity("Bernoulli.normal.rxprnSet", function(object)
{
	nam.ok <- sameNames(object, object@alternative, object@param)
	oks <- c(nam = nam.ok)
	ok <- all(oks==T)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})


setClass("XprnSet", representation("xprnSet"))
setValidity("XprnSet", function(object)
{
	all(as.numeric(exprs(object)) >= 0, na.rm = TRUE)
})
setMethod("print", signature(x = "XprnSet"), function(x)
{
  message(class(x), "; ratio or nonnegative intensity data:\n")
  print(as(x, "ExpressionSet"))
})
setMethod("exprs", signature(object = "XprnSet"), function(object)
{
  Matrix(exprs(as(object, "ExpressionSet")))
})
setMethod("logb", signature(x = "XprnSet", base = "missing"), function(x, base)
{
  logb(as(x, "xprnSet"))
})

coerceExpressionSet <- function(from, to.fun)
{
	assert.is(from, "ExpressionSet")
	assert.is(to.fun, "function")
	to.fun(phenoData = pData(from), exprs = exprs(from),featureData=fData(from))#April 2013: marta added featureData=fData(from)
}
setAs(from = "ExpressionSet", to = "xprnSet", function(from)
{
	coerceExpressionSet(from = from, to.fun = xprnSet)
})
setAs(from = "ExpressionSet", to = "XprnSet", function(from)
{
	coerceExpressionSet(from = from, to.fun = XprnSet) # as(as(from, "xprnSet"), "XprnSet")
})

setClass("xprnSetPair", representation(x = "xprnSet", y = "xprnSet"))
setValidity("xprnSetPair", function(object)
{
	cla.ok <- !is(object@x, "XprnSet") && !is(object@y, "XprnSet")
	fea.x <- featureNames(object@x)
	fea.y <- featureNames(object@y)
	fea.ok <- length(fea.x) == length(fea.x) && is.character(fea.x) && all(fea.x == fea.y)
	ok <- cla.ok && fea.ok
	if(!ok){printInvalid(object); browser()}
	ok
})
setAs(from = "xprnSetPair", to = "xprnSet", function(from)
{
	phenoData <- rbind(pData(from@x), pData(from@y))
	exprs <- cbind(exprs(from@x), exprs(from@y))
	xprnSet(phenoData = phenoData, exprs = exprs)
})

setMethod("featureNames", signature(object = "xprnSetPair"), function(object) # new 24 April 2008; for bias.r
{
  featureNames(object@x)
})
setMethod("names", signature(x = "xprnSetPair"), function(x) # new 24 April 2008; for bias.r
{
  featureNames(x)
})
setMethod("logb", signature(x = "xprnSetPair", base = "missing"), function(x, base)
{
	if(is(x@x, "XprnSet"))
	{
		x@x <- logb(x@x)
		x@y <- logb(x@y)
	}
	else
		warning("no log taken")
  x
})
setMethod("[", signature(x = "xprnSetPair", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@x <- x@x[i, ]
	x@y <- x@y[i, ]
	stopifnot(validObject(x))
	x
})

setClassUnion(name = "xprnSetObject", members = c("xprnSet", "xprnSetPair")) # used in estimate.s, biasEstimate
# setClassUnion("xprnObject", "xprnSetObject") # c("xprnSet", "xprnSetPair"))


setClass("xprnSetObjects", representation("list")) # new 28 April 2008
setValidity("xprnSetObjects", function(object)
{
	all(sapply(object, is, class2 = "xprnSetObject"))
})

setClass("xprnSetObjectPair", representation(training = "xprnSetObject", test = "xprnSetObject")) # new 24 April 2008; for bias.r
setValidity("xprnSetObjectPair", function(object)
{
	cla.ok <- TRUE
	fea.training <- featureNames(object@training)
	fea.test <- featureNames(object@test)
	fea.ok <- length(fea.training) == length(fea.training) && is.character(fea.training) && all(fea.training == fea.test)
	ok <- cla.ok && fea.ok
	if(!ok){printInvalid(object); browser()}
	ok
})
setMethod("featureNames", signature(object = "xprnSetObjectPair"), function(object) # new 24 April 2008; for bias.r
{
  featureNames(object@test)
})

setClass("smoothSpline", "list")
setValidity("smoothSpline", function(object)
{
	le.ok <- length(object) == 1
	if(!le.ok)
	{ printInvalid(object); browser()}
	ok <- le.ok && names(object) == "s3"
	if(is.na(ok) || !ok) warning("names lost")
	le.ok
})
setMethod("lines", signature(x = "smoothSpline"), function(x, ...)
{
  Lines(x = x, ...)
})
setMethod("plot", signature(x = "smoothSpline", y = "missing"), function(x, y, ...)
{
  Plot(x = x, ...)
})
setMethod("plot", signature(x = "smoothSpline", y = "smoothSpline"), function(x, y, ...)
{
  stopifnot(all())
  indices <- sort(union(s3(x, "x"), s3(y, "x")))
  get.y <- function(ss)
  {
#    ifelse(indices %in% s3(ss, "x"), s3(ss, "y"), as.numeric(NA))
    sapply(indices, function(i)
    {
    	s3(ss, "y")[which(s3(ss, "x") == i)]
    })
  }
  y.from.x <- get.y(x)
  y.from.y <- get.y(y)
  stopifnot(length(y.from.x) == length(y.from.y))
  plot(x = y.from.x, y = y.from.y, ...)
  abline(a = 0, b = 1, col = "blue")
})

setClass("Lowess", "smoothSpline")
setAs(from = "smoothSpline", to = "Lowess", function(from)
{
	lis <- as(from, "list")
	names(lis) <- "s3" # names(from)
	lo <- new("Lowess", lis)
	names(lo) <- "s3" # names(from)
	lo
})
setClass("movingAverage", "smoothSpline")
setAs(from = "smoothSpline", to = "movingAverage", function(from)
{
	movingAverage(x = from)
})
smoothData.classes <- c("Lowess", "movingAverage", "smoothSpline")
setClassUnion("smoothData", smoothData.classes)

setClass("movingLocation", representation(x = "numeric", y = "numeric"))
setMethod("lines", signature(x = "movingLocation"), function(x, ...)
{
  Lines(x = x, ...)
})
setMethod("plot", signature(x = "movingLocation", y = "missing"), function(x, y, ...)
{
  Plot(x = x, ...)
})

setClass("oneLeftOut", representation(training = "list", test = "list", noneLeftOut = "list"))
setValidity("oneLeftOut", function(object)
{
  len.ok <- length(object@training) >= 2 && length(object@training) == length(object@test) && length(object@noneLeftOut) == 1
  lens.ok <- all(length(noneLeftOut(object)) == sapply(object@test, length)) && all(length(noneLeftOut(object)) == sapply(object@training, length))
  cla.ok <- all(class(noneLeftOut(object)) == sapply(object@training, class) && class(noneLeftOut(object)) == sapply(object@test, class))
  ok <- len.ok && lens.ok && cla.ok
  if(!ok){message("invalid oneLeftOut"); browser()}
  ok
})
#XXX|:"xprnScore.p" and "xprnScore" are  not defined yet
#setMethod("sort", signature(x = "oneLeftOut", decreasing = "ANY"), function(x, decreasing = logical(0), verbose, ...)
#{
#	if(missing(verbose)) verbose <- FALSE
#	if(missing(decreasing) || length(decreasing) == 0) decreasing <- FALSE
#	nlo <- noneLeftOut(x)
#	if(is(nlo, "xprnScore.p"))#XXX|:"xprnScore.p" and "xprnScore", "P0" and "Fdr" are clases not defined yet
#		message('"sort", signature(x = "oneLeftOut", decreasing = "ANY") not yet implemented for ', class(nlo))
#	else {
#		sort(nlo,decreasing=decreasing)
#	}
#	#{
#	#	#sorted.olo <- sort(oneLeftOut(object = x, FUN = function(score.p)
#	#	#{
#	#	#	as.xprnScore(from = score.p, caller = '"sort", signature(x = "oneLeftOut", decreasing = "ANY")')
#	#	#}))
#	#	#sorted.olo.ok <- is(sorted.olo, "oneLeftOut") && is(noneLeftOut(sorted.olo), "xprnScore")
#	#	#if(!sorted.olo.ok)
#	#	#{ message("bad sorted.olo"); browser()}
#	#	#oneLeftOut(object = x, anotherLeftOut = sorted.olo, FUN = function(score.p, sorted.score)
#	#	#{
#	#	#	args.ok <- is(score.p, "xprnScore.p") && is(sorted.score, "xprnScore")
#	#	#	if(!args.ok)
#	#	#	{ message("list element of wrong class"); browser()}
#	#	#	unsorted.p <- score.p@p
#	#	#	if(is.null(names(unsorted.p)) || !sameNames(unsorted.p, sorted.score, order.sensitive = FALSE))
#	#	#	{ message("name inconsistency A"); browser()}
#	#	#	sorted.p <- unsorted.p[names(sorted.score)]
#	#	#	sorted.p@sorted <- TRUE
#	#	#	if(!sameNames(sorted.score, sorted.p, order.sensitive = TRUE))
#	#	#	{ message("name inconsistency B"); browser()}
#	#	#	xprnScore.p(sorted.score, p = sorted.p)
#	#	#})
#	#}
#	#else if(is(nlo, "xprnScore") || is(nlo, "estimate") || is(nlo, "sortableEstimate")) # classes defined in reprod.r and estimate.r
#	#{
#	#	get.ord <- function(object)
#	#	{
#	#		Or <- try(Order(object = object, decreasing = decreasing, ...))
#	#		if(is(Or, "try-error"))
#	#		{ message("bad Or"); browser()}
#	#		Or
#	#	}
#	#  for(i in 1:length(x@test))
#	#  {
#	#    trainingSco <- x@training[[i]]
#	#    testSco <- x@test[[i]]
#	#    stopifnot(is(trainingSco, class(nlo)) && is(testSco, class(nlo)) && length(trainingSco) == length(testSco))
#	#  	ord <- get.ord(object = trainingSco)
#	#  	x@training[[i]] <- trainingSco[ord]
#	#  	x@test[[i]] <- testSco[ord]
#	#  	stopifnot(all(names(x@training[[i]]) == names(x@test[[i]])))
#	#  }
#	#  nloSco <- noneLeftOut(x)
#	#  ord.nlo <- get.ord(object = nloSco)
#	#  x@noneLeftOut <- list(nloSco[ord.nlo])
#	#  if(all(names(nloSco) == names(noneLeftOut(x))))
#	#  { message("sort failed to change order of names"); browser()}
#	#  else if(verbose)
#	#		print(data.frame(names(nloSco), names(noneLeftOut(x)))[1:10, ])
#	#  x
#	#}
##	else
###	{
###
###	}
#	#{
#	message('"sort", signature(x = "oneLeftOut", decreasing = "ANY") not yet implemented for ', class(nlo)); browser()
#		#}
#})
#setMethod("plot", signature(x = "oneLeftOut", y = "missing"), function(x, y, call.sort, main, ...)
#{
#  nlo <- noneLeftOut(x)
#  if(missing(call.sort))
#  {
#    call.sort <- is(nlo, "xprnScore")
#    if(!is.factor(nlo))
#      message("reproducibility call.sort = ", call.sort)
#  }
#  if(missing(main)) main <- paste("n = ", length(x@test), "; based on ", class(nlo), sep = "")
#	if(is.factor(nlo) || is(nlo, "xprnScore"))
#	{
#		reprod <- reproducibility(x, call.sort = call.sort)
#		plot(x = reprod, main = main, ...)
#	}
#	else
#	  stop(paste('"plot", signature(x = "oneLeftOut", y = "missing") not yet implemented for ', class(noneLeftOut(x)), sep = ""))
#})

setClassUnion("Vector", c("numeric", "character", "logical"))
setClass("predictionError", representation("matrix", parameter = "Vector", parameter.lab = "character")) # parameter was "numeric" before 8 October 2007
setValidity("predictionError", function(object)
{
	len.ok <- ncol(object) == length(object@parameter) && length(object@parameter.lab) <= 1
	pe.ok <- !is.na(len.ok) && len.ok
	if(!pe.ok)
	{ printInvalid(object); browser()}
	pe.ok
})
setMethod("[", signature(x = "predictionError", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	rn <- rownames(x)
	x@.Data <- as(x, "matrix")[i, ]
	if(is.character(rn))
		rownames(x) <- rn[i]
	stopifnot(validObject(x))
	x
})
setMethod("plot", signature(x = "predictionError", y = "missing"), function(x, y, cumulative, indices, nfeatures, use.lower.indices, ylim, col, xlab, log, FUN, f, smooth, verbose, legend.x, ...)
{
	if(missing(legend.x)) legend.x <- "top"
	if(missing(verbose)) verbose <- FALSE
	if(verbose)
	{
		message("plotting ", class(x), " with cumulative ", if(missing(cumulative)) "[missing]" else cumulative, " and this for ...:")
		print(names(list(...)))
	}
	if(missing(log))
		log <- "y"
	if(missing(xlab))
		xlab <- "index or rank"
	if(missing(cumulative))
	{
		cumulative <- !missing(use.lower.indices)
		message("predictionError cumulative = ", cumulative)
	}
	if(missing(indices))
	{
		if(missing(nfeatures) || length(nfeatures) == 0)
		{
			nfeatures <- if(cumulative) min(2000, nrow(x)) else nrow(x)
			message("nfeatures = ", nfeatures)
		}
		if(missing(use.lower.indices))
		{
			use.lower.indices <- if(cumulative) TRUE else logical(0)
			message("use.lower.indices = ", use.lower.indices)
		}
		nf <- ceiling(min(nfeatures, if(length(use.lower.indices) == 0) nrow(x) / 2 else nrow(x)))
		lower <- 1:nf
		upper <- (nrow(x) + 1 - nf):nrow(x)
		indices <- if(length(use.lower.indices) == 0)
			union(lower, upper)
		else if(use.lower.indices)
			lower
		else
			upper
	}
	if(missing(smooth))
	{
		smooth <- if(cumulative)
			FALSE
		else
		{
			if(missing(f))
			{
				f <- 1/100
				if(missing(smooth) || !identical(smooth, FALSE))
					message("f = ", f)
			}
			stopifnot(is.numeric(f) && length(f) == 1)
			smoothDataFUN(f = f)
		}
	}
	smooth.x <- function(vec)
	{
		stopifnot(is.numeric(vec))
		if(is.logical(smooth) && !smooth)
			vec # cf. plot biasEstimates
		else
			smooth(vec)
	}
	x <- smooth.x(x)
	if(missing(FUN))
	{
		FUN.name <- if(cumulative) "identity" else "exp_root"
		message(class(x), " FUN = ", FUN.name)
		FUN <- eval(parse(text = FUN.name))
	}
	stopifnot(is.function(FUN))
	sub.y.mat <- if(cumulative)
	{
		if(verbose) message("calling accumulate with indices from ", min(indices), " to ", max(indices), ".")
		accumulate(x, indices = indices)
	}
	else
		x[indices, ]
	if(missing(ylim))
		ylim <- FUN(range(as.numeric(sub.y.mat), na.rm = TRUE))
	graph <- function(fun, indices, sub.y, ...)
	{
		stopifnot(length(sub.y) <= nrow(x))
#		sub.y <- get.sub.y(indices = indices, y = y)
		stopifnot(length(sub.y) == length(indices))
		fun(x = indices, y = FUN(sub.y), ...)
	}
	pch <- 1:length(x@parameter)
	if(missing(col))
		col <- pch
	if(length(col) == 1)
		col <- rep(col, length(x@parameter))
	stopifnot(length(col) == length(pch))
	arglis <- list(...)
	for(j in 1:length(x@parameter))
	{
		param <- x@parameter[j]
		fun <- if(j == 1)
			function(...)
			{
				do.call("plot", c(arglis, list(xlab = xlab, ylab = "prediction error", log = log, ylim = ylim, ...)))
			}
		else
			points
		graph(fun = fun, indices = indices, sub.y = sub.y.mat[, j], col = col[j], pch = pch[j])
	}
	legend(legend = if(length(x@parameter.lab) == 1) paste(x@parameter.lab, x@parameter, sep = " = ") else x@parameter, x = legend.x, col = col, pch = pch)
})

# setClass("general.numeric")
# setValidity("general.numeric", function(object)
# {
# 	ok <- is(object, "numeric")
# 	if(!ok)
# 	{ printInvalid(object); browser()}
# 	ok
# })

if(!isClass("numericEstimate"))
	setClassUnion("numericEstimate", c("numeric")) # modified in the file "estimate.r"

setClass("biasEstimate", representation("numeric", uncorrected = "numericEstimate", sorted = "logical", jackknife = "logical", Rank = "Numeric", Weight = "Scalar")) # changed 16 April 2008
setValidity("biasEstimate", function(object)
{
	ra <- object@Rank
	ra.ok <- length(ra) == 0 || (object@sorted && length(ra) == length(object) && all(names(ra) == names(object)) && all(ra >= 1))
	cla.ok <- is(object@uncorrected, "numeric")
	is.boo.ok <- function(boo){length(boo) == 1}
	boo.ok <- is.boo.ok(object@sorted) && is.boo.ok(object@jackknife)
	len.ok <- length(object) == length(object@uncorrected)
	nam.ok <- all(names(object) == names(object@uncorrected))
	ok <- ra.ok && cla.ok && boo.ok && !is.na(nam.ok) && nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("names", signature(x = "biasEstimate"), function(x){names(x@uncorrected)})
setMethod("annotation", signature(object = "biasEstimate"), function(object)
{
	ann <- try(annotation(object@uncorrected), silent = TRUE)
	if(is(ann, "try-error"))
		"uncorrected estimate"
	else
		ann
})
setAs(from = "biasEstimate", to = "numeric", function(from)
{
	vec <- as.numeric(from)
	names(vec) <- names(from)
	vec
})
setMethod("[", signature(x = "biasEstimate", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  bias <- as(x, "numeric")[i]
  uncorrected <- x@uncorrected[i]
  Ra <- if(x@sorted)
  	Rank(x)[i]
  else
  	x@Rank # Numeric(numeric(0))
  new.biasEstimate(bias = bias, uncorrected = uncorrected, sorted = x@sorted, jackknife = x@jackknife, Rank = Ra)
})
setMethod("plot", signature(x = "biasEstimate", y = "missing"), function(x, y, name, sub, xlab, ...)
{
	if(missing(name))
	{
		get.type <- function() # requires fdr.s code
		{
			x@uncorrected@estimator@type
		}
		type <- try(get.type(), silent = TRUE)
		name <- if(is(type, "try-error"))
			"uncorrected"
		else
			type
		message("name = ", name)
	}
	new.x <- if(name == "pvalue.z")
		qvalue(x)
	else if(name %in% slotNames(x))
		slot(x, name = name)
	else
		stop('bad name in "plot", signature(x = "biasEstimate", y = "missing")')
	if(missing(sub))
	{
		sub <- if(name == "pvalue.z")
			""
		else
			name
	}
	if(missing(xlab))
	{
		xlab <- if(name == "pvalue.z")
			"q-value"
		else
			annotation(x)
	}
	plot(x = new.x, y = x, xlab = xlab, ylab = "estimate", sub = sub, ...)
})
setMethod("plot", signature(x = "numeric", y = "biasEstimate"), function(x, y, xlab, ylab, smooth, call.browser, save.space, f, include.estimate, main, call.par, include.rough, ...)
{
	if(!y@sorted)
	{
		sort.y <- function()
		{
			message("sorting ", class(y), " on ", date())
			sorted.y <- sort(y)
			if(!sorted.y@sorted)
			{ message("bad sorted.y"); browser()}
			sorted.y
		}
		y <- sort.y()
	}
	if(missing(include.estimate))
	{
		include.estimate <- TRUE
		message("include.estimate = ", include.estimate)
	}
	if(missing(include.rough))
	{
		include.rough <- !include.estimate
		message("include.rough = ", include.rough)
	}
	if(missing(save.space))
	{
		save.space <- TRUE
		message("save.space = ", save.space)
	}
	stopifnot(length(x) == length(y) && is.character(names(x)) && is.character(names(y)))
	if(missing(call.par))
		call.par <- TRUE
	same.names <- function()
	{
		length(x) == length(y) && is.character(names(x)) && is.character(names(y)) && all(names(x) == names(y))
	}
	if(!same.names())
		x <- x[names(y)]
	stopifnot(same.names())
	if(missing(call.browser)) call.browser <- FALSE
	if(call.browser) browser()
	if(missing(smooth))
		smooth <- NULL
	if(missing(f))
		f <- NULL
	if(include.rough)
		y.rough <- y
	y <- Smooth(object = y, smooth = smooth, f = f)
	if(include.estimate && call.par) par(mfrow = c(2, 2))
	uncorrected.pch <- 1
	corrected.pch <- 2
	uncorrected.col <- "orange"
	corrected.col <- "blue"
	if(missing(main))
		main <- annotation(y)
	graph <- function(FUN, y.vec, ...)
	{
		FUN(x = as.numeric(x), y = as.numeric(y.vec), ...)
	}
	if(missing(ylab))
		ylab <- "estimate"
	if(include.estimate)
	{
		graph(FUN = plot, y.vec = y@uncorrected, xlab = xlab, ylab = ylab, pch = uncorrected.pch, col = uncorrected.col, main = main)
		graph(FUN = points, y.vec = corrected(y), pch = corrected.pch, col = corrected.col)
		if(!save.space)
			blank.plot()
		legend(x = if(save.space) "right" else "left", legend = c("uncorrected estimate", "corrected estimate"), pch = c(uncorrected.pch, corrected.pch), col = c(uncorrected.col, corrected.col), bg = "white")
	}
	rough.col <- "black"
	graph(FUN = plot, y.vec = if(include.rough) y.rough else y, xlab = xlab, ylab = "estimated bias", main = main, col = if(include.rough) rough.col else "black", ...)
	if(include.rough)
	{
		if(all(y.rough == y, na.rm = TRUE))
			warning("unsmoothed points do not differ from smoothed points")
		graph(FUN = points, y.vec = y, type = "l", col = "orange", lwd = 2)
	}
})
setMethod("plot", signature(x = "biasEstimate", y = "biasEstimate"), function(x, y, call.par, ...)
{
	if(missing(call.par))
		call.par <- TRUE
	if(call.par)
		par(mfrow = c(2, 2))
	xlab <- annotation(x)
	ylab <- annotation(y)
	graph <- function(FUN, main)
	{
		get.vec <- function(object){as.numeric(FUN(object))}
		plot(x = get.vec(x), y = get.vec(y), xlab = xlab, ylab = ylab, main = main, ...)
		abline(a = 0, b = 1, col = "blue")
	}
	graph(FUN = function(object){object@uncorrected}, main = "uncorrected estimates")
	graph(FUN = corrected, main = "corrected estimates")
	graph(FUN = identity, main = "bias")
})

setClass("Density", representation(s3 = "list", annotation = "character"))
setValidity("Density", function(object)
{
	s3.ok <- length(object@s3) == 1 && names(object@s3) == "s3"
	ann.ok <- length(object@annotation %in% c(0, 1))
	ok <- s3.ok && ann.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "Density", y = "missing"), function(x, y, main, ...)
{
	if(missing(main)) main <- annotation(x)
	plot(s3(x), main = main, ...)
})
setAs(from = "Density", to = "function", function(from)
{
	den <- s3(from)
	dom <- domain(from)
	function(v)
	{
		adens <- approxfun(x = den$x, y = den$y)(v)
		stopifnot(length(v) == length(adens))
		dens <- try(ifelse(v >= dom[1] & v <= dom[2], adens, 0))
		if(is(dens, "try-error"))
		{ message("bad dens"); browser()}
		dens
	}
})

setClass("functions", representation("list"))
setValidity("functions", function(object)
{
	ok <- length(object) == 0 || all(sapply(object, is, class2 = "function"))
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})

setMethod("plot", signature(x = "list", y = "missing"), function(x, y, call.par, ...) # based on vplot3
{
	if(missing(call.par))
		call.par <- length(x) > 1
	if(call.par)
		par(mfrow=c(2,2))
	for(i in 1:length (x))
	{
		message("plotting ", class (x[[i]]))
		plot(x=x[[i]], ...)
	}
})

setAs(from = "ANY", to = "list", function(from)
{
	nam <- slotNames(from)
	lis <- lapply(nam, function(name){slot(from, name = name)})
	names(lis) <- nam
	lis
})

setClass("binomMetalevel", representation(valid = "Numeric", nonconservative = "Numeric", corrected = "Numeric", x = "Numeric", size = "Numeric", prob1 = "Numeric", prob2 = "Numeric"))
setValidity("binomMetalevel", function(object)
{
	probs.ok <- length(object@prob1) == length(object@prob2) && all(object@prob2 >= object@prob1)
	lis <- list(object@valid, object@nonconservative, object@corrected)
	lens <- sapply(lis, length)
	len.ok <- all(length(object) == lens)
	len1.ok <- length(object@prob1) == 1 || length(object@x) == 1
	na.ok <- !any(is.na(c(object@valid, object@nonconservative)))
	oks <- c(probs = probs.ok, len = len.ok, len1 = len1.ok, na = na.ok)
	ok <- all(oks==T)
	if(!ok)
	{ printInvalid(object); print(oks); browser()}
	ok
})
setAs(from = "binomMetalevel", to = "numeric", function(from) # return(maxEntropy(from)) prior to 091026
{
	meanOverConvexSet(from)
})
setMethod("length", signature(x = "binomMetalevel"), function(x){max(length(x@x), length(x@prob1))})
setMethod("print", signature(x = "binomMetalevel"), function(x, ...)
{
	print(data.frame(valid = as.numeric(x@valid), nonconservative = as.numeric(x@nonconservative), x = as.numeric(x@x), size = as.numeric(x@size)))
})
setMethod("plot", signature(x = "binomMetalevel", y = "missing"), function(x, y, call.maxEntropy = FALSE, ...)
{
	vs.size <- length(x@prob1) == 1
	X <- if(vs.size)
		x@size
	else
	{
		stopifnot(length(x@prob1) == length(x) && length(x) >= 1)
		x@prob2 - x@prob1
	}
	graph <- function(fun, y, ...){fun(x = X, y = y, ...)}
	pch.valid <- ">"
	pch.nc <- "<"
	pch.c <- "h"
	pch.maxent <- "\\"
	col.valid <- col.nc <- col.c <- "black"
	col.maxent <- "gray"
	xlab <- if(vs.size)
		"number of trials"
	else
		"hypothesis interval width for binomial parameter"
	ylab <- if(vs.size)
		paste("confidence that ", x@prob1, " < p < ", x@prob2, sep = "")
	else
		"confidence that parameter is in hypothesis interval"
	graph(fun = plot, y = x@valid, xlab = xlab, ylab = ylab, ylim = c(0, 1), pch = pch.valid, log = if(vs.size) "x" else "", ...)
	graph(fun = points, y = x@nonconservative, pch = pch.nc)
	graph(fun = points, y = x@corrected, pch = pch.c)
	graph(fun = lines, y = if(call.maxEntropy) maxEntropy(x) else as(x, "numeric"), col = col.maxent) # , pch = pch.maxent
	legend.maxent <- if(call.maxEntropy) "maximum entropy" else "mean over convex set"
	legend(legend = c("based on nonconservative CI", "based on valid CI", "based on half-correction", legend.maxent), pch = c(pch.nc, pch.valid, pch.c, pch.maxent), x = if(vs.size) "bottomright" else "topleft", col = c(col.nc, col.valid, col.c, col.maxent))
	if(vs.size)
		abline(h = 1, lty = "dashed")
})

setClass("argmax", representation("scalar", approx.max = "scalar", from = "scalar", to = "scalar"))
setClass("invert", representation("argmax", value = "scalar", approx.value = "scalar", fun = "function"))

setClass("ttest", representation(pvalue = "Numeric", stat = "numeric", df = "Numeric", alternative = "character", level1 = "character", level2 = "character"))
setValidity("ttest", function(object)
{
	len.ok <- all(length(object) == c(length(object@pvalue), length(object@df)))
	len1.ok <- all(1 == sapply(list(object@alternative, object@level1, object@level2), length))
	nam.ok <- sameNames(object@pvalue, object@stat, object@df)
	ok <- len.ok && len1.ok && nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("length", signature(x = "ttest"), function(x)
{
	length(x@stat)
})
setMethod("names", signature(x = "ttest"), function(x)
{
	names(x@stat)
})
setMethod("[", signature(x = "ttest", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@pvalue <- x@pvalue[i]
  x@stat <- x@stat[i]
  x@df <- x@df[i]
  x
})

setClass("rParameter", representation("functions", ann = "character"))

setClass("nonparametricBootstrap", representation("matrix", original = "numeric"))
setValidity("nonparametricBootstrap", function(object)
{
	len.ok <- length(object@original) %in% c(0, nrow(object))
	nam.ok <- length(object@original) == 0 || all(rownames(object) == names(object@original))
	ok <- len.ok && nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("mean", signature(x = "nonparametricBootstrap"), function(x, na.rm = TRUE)
{
	m <- sapply(1:nrow(x), function(i){mean(x[i, ], na.rm = na.rm)})
	if(nrow(x) != length(m))
	{ message("nrow(x) != length(m)"); browser()}
	names(m) <- rownames(x)
	m
})

setClass("CI", representation("matrix", expectation = "numeric", zeroCI = "numeric", point = "numeric"))
setValidity("CI", function(object)
{
	point.ok <- length(object@point) == 0 || (length(object@point) == nrow(object) && all(rownames(object) == names(object@point)))
	len.ok <- nrow(object) == length(object@expectation) && ncol(object) == 2
	nam.ok <- all(rownames(object) == names(object@expectation))
	ord.ok <- all(object[, 1] <= object[, 2])
	zeroCI.ok <- try(nrow(object) == length(object@zeroCI) && all(rownames(object) == names(object@zeroCI))) # backward-compatible
	oks <- c(point = point.ok, len = len.ok, nam = nam.ok, ord = ord.ok, zero = (is(zeroCI.ok, "try-error") || zeroCI.ok))
	ok <- all(oks==T)
	if(is.na(ok) || !ok)
	{ printInvalid(object); print(oks); browser()}
	ok
})
setMethod("names", signature(x = "CI"), function(x)
{
	rownames(x)
})
setAs(from = "list", to = "CI", function(from)
{
	assert.are(from, "CI")
	get.limit <- function(j)
	{
		get.scalar <- function(ci)
		{
			sca1 <- ci[, j]
			stopifnot(length(sca1) == 1)
			sca1
		}
		sapply(from, get.scalar)
	}
	object <- cbind(get.limit(j = 1), get.limit(j = 2))
	get.slot <- function(name)
	{
		get.scalar <- function(ci)
		{
			sca <- slot(ci, name = name)
			if(length(sca) > 1)
			{ message("length(sca) > 1"); print(sca); browser()}
			sca
		}
		vec <- sapply(from, get.scalar)
		if(is.numeric(vec))
			vec
		else
			numeric(0)
	}
	expectation <- get.slot(name = "expectation")
	zeroCI <- get.slot(name = "zeroCI")
	point <- get.slot(name = "point")
	new.CI(object = object, expectation = expectation, zeroCI = zeroCI, point = point)
})
setAs(from = "CI", to = "data.frame", function(from)
{
  mat <- cbind(as(from, "matrix"), expectation = from@expectation, point = from@point)
  as.data.frame(mat)
})
#setMethod("[", signature(x = "CI", i = "ANY", j = "logical", drop = "missing"), function(x, i, j, drop)
#{
#  new.CI(object = as(x, "matrix")[i, ], expectation = x@expectation[i], point = x@point[i])
#	mat <- x@.Data[i, ]
#	if(!is(mat, "matrix"))
#	{ message("mat is ", class(mat)); browser()}
#	x@.Data <- mat
#	x@expectation <- x@expectation[i]
#	x@point <- x@point[i]
#	stopifnot(validObject(x))
#	x
#})
#setMethod("subset", signature(x = "CI", subset = "ANY", select = "missing", drop = "missing"), function(x, subset, select, drop)
#{
#  new.CI(object = as(x, "matrix")[i, ], expectation = x@expectation[i], point = x@point[i])
#	mat <- x@.Data[subset, ]
#	if(!is(mat, "matrix"))
#	{ message("* mat is ", class(mat)); browser()}
#	x@.Data <- mat[subset, ]
#	x@expectation <- x@expectation[subset]
#	x@point <- x@point[subset]
#	stopifnot(validObject(x))
#	x
#})

setClass("pointEstimate", representation(expectation = "numeric", zeroCI = "numeric", point = "numeric"))
setValidity("pointEstimate", function(object)
{
	point.ok <- length(object@point) == 0
	zeroCI.ok <- try(length(object@expectation) == length(object@zeroCI) && all(names(object@expectation) == names(object@zeroCI)))
	oks <- c(point = point.ok, zero = (is(zeroCI.ok, "try-error") || zeroCI.ok))
	ok <- all(oks==T)
	if(is.na(ok) || !ok)
	{ printInvalid(object); print(oks); browser()}
	ok
})
setMethod("length", signature(x = "pointEstimate"), function(x)
{
	length(x@zeroCI)
})
setMethod("names", signature(x = "pointEstimate"), function(x)
{
	names(x@expectation)
})
setAs(from = "list", to = "pointEstimate", function(from)
{
	assert.are(from, "pointEstimate")
	get.slot <- function(name)
	{
		get.scalar <- function(ci)
		{
			sca <- slot(ci, name = name)
			if(length(sca) > 1)
			{ message("length(sca) > 1"); print(sca); browser()}
			sca
		}
		vec <- sapply(from, get.scalar)
		if(is.numeric(vec))
			vec
		else
			numeric(0)
	}
	expectation <- get.slot(name = "expectation")
	zeroCI <- get.slot(name = "zeroCI")
	point <- get.slot(name = "point")
	new.pointEstimate(expectation = expectation, zeroCI = zeroCI, point = point)
})

setClassUnion(name = "posteriorEstimate", members = c("CI", "pointEstimate"))
setMethod("plot", signature(x = "posteriorEstimate", y = "missing"), function(x, y, add = FALSE, i = 1:length(x), col = "black", name, ...)
{
	if(missing(name))
	{
		name <- if(is.bootstrap(x))
			"point"
		else
			"zeroCI"
	}
	point <- slot(x, name = name)[i]
	expectation <- x@expectation[i]
	graph <- function(fun, ...)
	{
		fun(x = point, y = expectation, col = col, ...)
	}
	pl <- try(if(add)
		graph(fun = points)
	else
		graph(fun = plot, xlab = "point estimate", ylab = "expectation", ...))
	if(is(pl, "try-error"))
	{ message("plot err"); print(c(length(point), length(expectation))); browser()}
	abline(a = 0, b = 1, col = "gray")
})
setMethod("plot", signature(x = "posteriorEstimate", y = "ttest"), function(x, y, neg.col = "gray", pos.col = "black", ...)
{
	y <- y[names(x)]
	pos.boo <- y@stat > 0
	graph <- function(i, add, col, ...)
	{
		plot(x = x, add = add, i = i, col = col, ...)
	}
	graph(i = !pos.boo, add = FALSE, col = neg.col, xlim = Lim(x@point), ylim = Lim(x@expectation), ...)
	graph(i = pos.boo, add = TRUE, col = pos.col)
	legend(legend = c("left tail", "right tail"), pch = 1, x = "bottomright", col = c(neg.col, pos.col))
})

setClass("Coverage", representation("numeric", interest.value = "numeric", covered.value = "matrix", nominal.rate = "Scalar", base = "Scalar"))
setValidity("Coverage", function(object)
{
	ncol.ok <- ncol(object@covered.value) %in% c(0, 1, length(object))
	len.ok <- length(object@covered.value) == 0 || ncol(object@covered.value) == length(object@interest.value)
	rate.ok <- length(object) == length(object@interest.value)
	ok <- ncol.ok && len.ok && rate.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "Coverage", y = "missing"), function(x, y, ...)
{
	X <- x@interest.value
	Y <- percent(x)
	plot(X, Y, xlab = "parameter value", ylab = "coverage rate (%)", ...) # xlim = Lim(), ylim = Lim(), ...
})

setClass("trivialFunction", representation("function", min_param = "scalar", max_param = "scalar"))
setClassUnion("Function", "trivialFunction") # later extended to "CDF" and/or "dposterior"
setMethod("plot", signature(x = "Function", y = "missing"), function(x, y, xlim, ylim = if(is(x, "CDF")) c(0, 1) else NULL, h, save.time = default(FALSE, "save.time"), length.out, ...)
{
	if("ylim" %in% names(list(...))) { message("ylim in ... (", class(x), ", missing)"); browser()}
	if(missing(h))
	{
		h <- default(if(is(x, "pposterior")) numeric(0) else c(0.25, 0.5, 0.75), "h")
		message("plotting Function of class ", class(x), "\n")
	}
	if(missing(xlim))
		xlim <- domain(x)
	if(missing(length.out))
	{	length.out <- if(save.time) 100 else 1000}
	by <- (xlim[2] - xlim[1]) / length.out
	param.value <- seq(xlim[1] + by, xlim[2] - by, by = by)
	message("plotting Function over domain of ", length(param.value), " values on", date())
	Y <- plot(x = param.value, y = x, xlim = xlim, ylim = ylim, save.time = save.time, ...)
	hvline(x = param.value, y = Y, h = h, save.time = save.time)
})
# setMethod("nrow", signature(x = "Coverage
#XXX|:"xprnScore.p" and "xprnScore", "P0" and "Fdr" are clases not defined yet
# near end of file:

#source(file = "data.s") # functions
