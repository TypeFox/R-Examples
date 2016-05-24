# congruity.r created by David Bickel on 1 July 2009.
# indicatorRisk changed to probabilisticRisk on 090717 15:20

#source("cdf.r")
#library(locfdr)
#library(fBasics)

setClass("probMass", representation(x = "numeric", prob = "Numeric"))
setValidity("probMass", function(object)
{
	ok <- length(object@x) > 0 && length(object@x) == length(object@prob) && sum(object@prob) == 1
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setAs(from = "probMass", to = "function", function(from)
{
	function(){sample(x = from@x, prob = from@prob, replace = TRUE, size = 1)}
})

setClass("marginal.zz", representation("numeric", pre.marginal = "numeric")) # added 090721
setValidity("marginal.zz", function(object)
{
	ok <- length(object) == length(object@pre.marginal)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})

setClass("alt", representation("character"))
setValidity("alt", function(object)
{
	len.ok <- length(object) %in% c(0, 1)
	val.ok <- length(object) == 0 || object %in% c("less", "greater", "Greater", "two.sided")
	ok <- len.ok && val.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setClass("empiricalNull", representation(PDF = "function", PDF0 = "function", CDF0 = "CDF", p0 = "Scalar", s3 = "list"))
setValidity("empiricalNull", function(object)
{
	ok <- TRUE # object@p0 < 1.01 (cf. max.p0)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setClass("cvalue", representation(zz = "numeric", s3FUN = "function", arglis = "list"))
setValidity("cvalue", function(object)
{
	alternative <- alt(object)
	alt.ok <- is(alternative, "alt") && (length(alternative) == 0 || alternative != "two.sided")
	cla.ok <- !is(object@zz, "Numeric")
	ok <- alt.ok && cla.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("[", signature(x = "cvalue", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@zz <- x@zz[i]
	stopifnot(validObject(x))
	x
})
setMethod("plot", signature(x = "cvalue", y = "missing"), function(x, y, main = "c-values", add.col = NULL, call.hist = is.null(add.col), ...)
{
	if(call.hist)
	{
		par(mfrow = c(2, 2))
		hist(x@zz, main = main, xlab = "z-transform of certainty level")
	}
	con <- sort(congruity(x))
	graph <- function(fun, ...)
	{
		fun(x = con, y = (1:length(con)) / (length(con) + 1), ...)
	}
	pl <- try(graph(fun = plot, xlab = "certainty level", ylab = "empirical CDF", main = main, ...))
	if(is(pl, "try-error"))
	{ message("cannot plot certainty level"); browser()}
	if(!is.null(add.col))
		graph(fun = points, col = add.col, ...)
	abline(a = 0, b = 1, col = "gray")
})
setMethod("plot", signature(x = "cvalue", y = "cvalue"), function(x, y, main = "c-values", ...)
{
	get.vec <- function(object)
	{
		con <- congruity(object)
		if(is.null(names(con)))
		{
			if(length(x) == length(y))
				names(con) <- make.names(1:length(con))
			else
				stop("unnamed x and/or y of unequal length")
		}
		con
	}
	X <- get.vec(x)
	Y <- get.vec(y)
	get.names <- function(object)
	{
		boo <- is.finite(object)
		if(!any(boo))
		{ message("no c-values to plot"); browser()}
		names(object[boo])
	}
	nam <- intersect(get.names(X), get.names(Y))
	if(length(nam) == 0)
		stop("nothing to plot")
	plot(X[nam], Y[nam], main = main, ...)
})
setMethod("names", signature(x = "cvalue"), function(x)
{
	names(x@zz)
})
setMethod("length", signature(x = "cvalue"), function(x)
{
	length(x@zz)
})
setAs(from = "cvalue", to = "numeric", function(from)
{
	x <- try(from@zz)
	if(is(x, "try-error"))
	{ message("cannot coerce cvalue to numeric"); browser()}
	boo <- !is.na(x) # is.finite before 090711
	if(!all(boo))
	{
		message("  coersion from ", class(from), " retaining ", sum(boo), " present z values of ", length(x), " original z values on ", date())
		x <- x[boo]
	}
	x
})
setAs(from = "ttest", to = "cvalue", function(from)
{
	new.cvalue(pvalue = from@pvalue, s3FUN = t.test, arglis = list(alternative = from@alternative))
})
setMethod("Combine", signature(x = "cvalue", y = "cvalue"), function(x, y, xlab = stop("no xlab"), ylab = stop("no ylab"), ...) # like x = "ttest", y = "ttest"
{
	assert.is(xlab, "character")
	assert.is(ylab, "character")
	stopifnot(length(xlab) == 1)
	stopifnot(length(ylab) == 1)
	stopifnot(xlab != ylab)
	same <- function(name)
	{
		identical(slot(x, name), slot(y, name))
	}
	stopifnot(same("s3FUN"))
	stopifnot(same("arglis"))
	get.z <- function(object, lab)
	{
		z <- object@zz
		names(z) <- sapply(names(z), function(na){paste(lab, na, sep = ".")})
		z
	}
	p1 <- get.z(x, lab = xlab)
	p2 <- get.z(y, lab = ylab)
	zvalue <- c(p1, p2)
	dup <- duplicated(names(zvalue))
	if(any(dup))
	{ message("duplication in names(zvalue):"); print(dup); browser()}
	x@zz <- zvalue
	x
#	new.cvalue(pvalue = pvalue, s3FUN = same("s3FUN"), arglis = list(alternative = same("alternative")))
})


setClass("ran.cvalue", representation("cvalue", param.sign = "numeric", p0 = "Scalar", mean = "scalar", sd = "Scalar", df = "Scalar", mean0 = "scalar", sd0 = "Scalar", df0 = "Scalar"))
setValidity("ran.cvalue", function(object)
{
	len.ok <- length(object) == length(object@param.sign)
	sign.ok <- all(object@param.sign %in% c(-1, 1))
	ok <- len.ok && sign.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "ran.cvalue", y = "missing"), function(x, y, ...)
{
	hist(x@zz, xlab = "z", xlim = c(-3, 3)) # sapply(1:1000, function(i){random.value(x)}))
})

setClass("adjusted.cvalue", representation("cvalue", empirical.null = "empiricalNull"))
setValidity("adjusted.cvalue", function(object)
{
	len.ok <- length(object) == length(lfdr(object@empirical.null))
	ok <- len.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "cvalue", y = "adjusted.cvalue"), function(x, y, statistic, statistic.name, null.congruity.range, statistic.range = NULL, call.ecdf = FALSE, ...)
{
	nam <- intersect(names(x), names(y))
	get.con <- function(object){congruity(object)[nam]}
	con <- get.con(x)
	acon <- get.con(y)
	if(!missing(statistic))
	{
		statistic <- statistic[nam]
		if(missing(statistic.name))
			statistic.name <- statistic.name
		graph <- function(fun, i, call.abline, adjust.ylim = FALSE, ...)
		{
			Y <- (acon - con)[i]
			finite.Y <- Y[is.finite(Y)]
			ylim <- if(adjust.ylim)
				c(-max(abs(finite.Y)), max(abs(finite.Y)))
			else
				range(finite.Y)
			fun(x = exp(statistic)[i], y = Y, ylim = ylim, ...)
			abline(h = 0, col = "gray")
		}
		use.ncr <- if(missing(null.congruity.range) || is.null(null.congruity.range))
			FALSE
		else if(is.numeric(null.congruity.range) && length(null.congruity.range) == 2)
			TRUE
		else
			stop("bad args in plot cvalue, adjusted.cvalue")
		figure <- function(col)
		{
			graph(fun = plot, i = nam, call.abline = TRUE, xlab = if(statistic.name == "mean") "expression ratio estimate" else paste("exp", statistic.name), ylab = paste("adjustment to certainty of ratio < 1"), log = "x", col = col, adjust.ylim = TRUE, ...)
		}
		all.col <- "gray"
		figure(col = if(use.ncr) all.col else "black")
		if(use.ncr)
		{
			get.names <- function(con, negate, sign = NULL)
			{
				boo <- con < null.congruity.range[1] | con > null.congruity.range[2]
				if(negate)
					boo <- !boo
				na <- names(con[boo])
				if(missing(sign) || is.null(sign))
					na
				else
				{
					boo <- if(sign == -1)
						con > 0.5
					else if(sign == +1)
						con < 0.5
					else
						stop("a bad sign")
					statistic.nam <- names(con[boo])
					intersect(na, statistic.nam)
				}
			}
			get.demoted.i <- function(sign = NULL)
			{
				intersect(get.names(con, negate = FALSE, sign = sign), get.names(acon, negate = TRUE, sign = sign))
			}
			demoted.i <- get.demoted.i()
			promoted.i <- intersect(get.names(con, negate = TRUE), get.names(acon, negate = FALSE))
			message(length(demoted.i), " genes demoted and ", length(promoted.i), " genes promoted")
			demoted.col <- "black"
			promoted.col <- "orange"
			title(main = paste(length(demoted.i), "genes adjusted to insignificance"))
			get.legend <- function(legend, col)
			{
				legend(legend = legend, col = col, x = "bottomleft", pch = 1, bty = "n")
			}
			if(length(promoted.i) == 0)
			{
				under.col <- "orange"
				over.col <- "black"
				under.i <- get.demoted.i(sign = -1)
				over.i <- get.demoted.i(sign = +1)
				graph(fun = points, i = under.i, call.abline = FALSE, col = under.col, ...)
				graph(fun = points, i = over.i, call.abline = FALSE, col = over.col, ...)
#				get.legend(legend = c(paste(length(demoted.i), "demoted genes"), paste("all", length(nam), "genes")), col = c(demoted.col, all.col))
				legend.legend <- function(prefix, con.threshold, col)
				{
					paste(if(prefix == "under") "ratio < 1" else "ratio > 1", " no longer ", 100 * con.threshold, "% certain", sep = "")
				}
				under.legend <- legend.legend(prefix = "under", con.threshold = null.congruity.range[2])
				over.legend <- legend.legend(prefix = "over", con.threshold = 1 - null.congruity.range[1])
				get.legend(legend = c(under.legend, over.legend, paste("all", length(nam), "genes")), col = c(under.col, over.col, all.col))
			}
			else
			{
				graph(fun = points, i = demoted.i, call.abline = FALSE, col = demoted.col, ...)
				graph(fun = points, i = promoted.i, call.abline = FALSE, col = promoted.col, ...)
				get.legend(legend = c("demoted genes", "promoted genes", "all genes"), col = c(demoted.col, promoted.col, all.col))
			}
		}
	} # end if(!missing(statistic))
	if(call.ecdf)
	{
		stop("not yet implemented")
		ecdf.col <- "gray"
		null.col <- "black"
		cdf0 <- y@CDF0
		X <- seq(cdf0@min_param, cdf0@max_param, length.out = 1e3)
		plot(x = X, y = as(cdf0, "function"), xlab = "nominal certainty level", ylab = "cumulative probability estimate", main = "estimated certainty level CDFs", col = null.col, ...)
		abline(0, 1, col = "gray")
		abline(v = 0.5, h = 0.5, col = "gray")
		plot(x = X, y = ecdf(con), col = ecdf.col, add = TRUE, ...)
#		plot(as(x, "cvalue"), add.col = ecdf.col, call.hist = FALSE)
	}
	else # before 090717
		plot(x = con, y = acon, xlab = "nominal certainty level", ylab = "adjusted certainty level", main = "null CDF estimate", ...)
	if(!missing(statistic.range) && !is.null(statistic.range)) # statistic.range is like mu.range of odd.s
	{
		stop("statistic.range not yet implemented")
	}
})

setClass("unconditional.cvalue", representation("cvalue", nominal.zz = "numeric")) # added 090712
setValidity("unconditional.cvalue", function(object)
{
	len.ok <- length(object) == length(object@nominal.zz)
	ok <- len.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setAs(from = "unconditional.cvalue", to = "function", function(from) # converts z-score to z-score that is N(0,1) under the true global null
{
	from@s3FUN
})
plot_unconditional_cvalue <- function(nom.zz, true.zz, fun, ...)
{
	assert.are(list(nom.zz, true.zz), "numeric")
	plot(x = nom.zz, y = true.zz, xlab = "nominal certainty z-score", ylab = "marginal certainty z-score", ...)
	lines(x = nom.zz, y = sapply(nom.zz, fun), col = "blue")
	abline(0, 1, col = "gray")
	abline(v = 0, h = 0, col = "gray")
}
setMethod("plot", signature(x = "numeric", y = "unconditional.cvalue"), function(x, y, ...)
{
	fun <- as(y, "function")
	nom.zz <- x
	true.zz <- fun(nom.zz) # z-score that is N(0,1) under the true global null
	plot_unconditional_cvalue(nom.zz = nom.zz, true.zz = true.zz, fun = as(y, "function"), ...)
})
setMethod("plot", signature(x = "unconditional.cvalue", y = "missing"), function(x, y, ...)
{
	true.zz <- x@zz # z-score that is N(0,1) under the true global null
	nom.zz <- x@nominal.zz
	plot_unconditional_cvalue(nom.zz = nom.zz, true.zz = true.zz, fun = as(x, "function"), ...)
})

setClass("statistic.cvalue", representation("cvalue", statistic = "numeric", statistic.name = "character", statistic.fun = "function"))
setValidity("statistic.cvalue", function(object)
{
	len.ok <- length(object@statistic.name) == 1
	nam.ok <- sameNames(object, object@statistic)
	ok <- len.ok && nam.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("[", signature(x = "statistic.cvalue", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	new.statistic.cvalue(object = as(x, "cvalue")[i], statistic = x@statistic[i], statistic.name = x@statistic.name, statistic.fun = x@statistic.fun)
})
setMethod("plot", signature(x = "statistic.cvalue", y = "adjusted.cvalue"), function(x, y, ...)
{
	tr <- try(plot_statistic.cvalue(x = x, acval = y, ...))
	if(is(tr, "try-error"))
	{ message("bad plot"); browser()}
})
plot_statistic.cvalue <- function(x, y, null.congruity.range = NULL, transformed = default(TRUE, "transformed"), statistic.range = NULL, call.par = TRUE, xlab, call.hist = default(FALSE, "call.hist"), acval, ...)
{
	assert.is(x, "statistic.cvalue")
	stopifnot(missing(y))
	if(call.par)
		par(mfrow = c(2, 2))
	x22 <- function(){if(!call.par) nx11()}
	plot.last <- function(){message("plot.last was not assigned")}
	if(missing(statistic.range) || is.null(statistic.range)) # statistic.range is like mu.range of odd.s
	{
		center <- if(transformed) 0 else 0.5
		get.con <- function(object){if(transformed) object@zz else congruity(object)}
		statistic <- x@statistic
		statistic.name <- x@statistic.name
		if(missing(xlab))
			xlab <- if(statistic.name == "mean") "expression ratio estimate" else paste("exp", statistic.name)
		if(missing(acval))
			acval <- try(adjusted.cvalue(x, plot = 0))
		assert.is(acval, "adjusted.cvalue")
		x.name <- "nominal"
		acval.name <- "adjusted"
		if((missing(null.congruity.range) || is.null(null.congruity.range)) && !call.hist)
		{
			plot_statistic <- function(cval, ylab)
			{
				con <- get.con(cval)
				if(transformed) ylab <- paste("transformed", ylab)
				plot(x = exp(statistic[names(con)]), y = con, xlab = xlab, ylab = ylab, log = "x", ...)
				abline(h = center, col = "gray")
			}
			get.lab <- function(name){paste(name, "congruity")}
			x22()
			plot_statistic(cval = x, ylab = get.lab(x.name))
			plot_statistic(cval = acval, ylab = get.lab(acval.name))
		}
		else if(is.numeric(null.congruity.range) && length(null.congruity.range) == 2)
		{
			if(call.hist)
			{
				histogram <- function(cval, main)
				{
					con <- get.con(cval)
					statistic <- statistic[names(con)]
					boo <- con < null.congruity.range[1] | con > null.congruity.range[2]
					X <- statistic[boo]
					hi <- try(hist(X, main = main, xlab = xlab, ylab = paste("number of genes out of", sum(is.finite(X)))))
					if(is(hi, "try-error"))
					{ message("no hist"); browser()}
				}
				get.main <- function(name){paste(name, " c outside [", null.congruity.range[1], ", ", null.congruity.range[2], "]", sep = "")}
				histogram(cval = x, main = get.main(x.name))
				histogram(cval = acval, main = get.main(acval.name))
			}
			else # used to generate "c10.emf" for "significance testing.lyx" via "figs 2009-07-17 408 33.pdf"
			{
				graph <- function(fun, cval, ...)
				{
					con <- get.con(cval)
					statistic <- statistic[names(con)]
					boo <- con < null.congruity.range[1] | con > null.congruity.range[2]
					Y <- if(transformed)
						con
					else
						con[boo]
					Y <- Y[is.finite(Y)]
					X <- exp(statistic[names(Y)])
					nonhist.main <- if(transformed)
						paste("all", length(X), "genes with 0 < p < 1")
					else
						paste("c outside [", null.congruity.range[1], ", ", null.congruity.range[2], "]", sep = "")
					fun(x = X, y = Y, xlab = xlab, ylab = paste(if(transformed) "z-transformed" else "underexpression", "certainty of ratio < 1"), main = nonhist.main, log = "x", ...)
				}
				cval.col <- "gray"
				acval.col <- "black"
				x22()
				graph(fun = plot, cval = x, col = cval.col, ...)
				graph(fun = points, cval = acval, col = acval.col, ...)
				abline(h = if(transformed) qnorm(null.congruity.range) else null.congruity.range, col = "gray")
				legend(legend = c("nominal certainty level", "adjusted certainty level"), col = c(cval.col, acval.col), x = "bottomleft", pch = 1, bty = "n")
				plot.last <- function(){plot(as(x, "cvalue"), call.hist = FALSE)}
			}
		}
		else
			stop("bad args in plot")
		plot(x = as(x, "cvalue"), y = acval, statistic = statistic, statistic.name = statistic.name, null.congruity.range = null.congruity.range, ...)
	} # end if(missing(statistic.range) || is.null(statistic.range))
	else
	{
		stop("statistic.range not yet implemented")
	}
	plot.last()
}
setMethod("plot", signature(x = "statistic.cvalue", y = "missing"), function(x, y, simple = default(x@statistic.name == "sign decision", "simple"), ...)
{
	if(!is.null(simple) && !simple)
	{
		tr <- try(plot_statistic.cvalue(x = x, ...))
		if(is(tr, "try-error"))
		{
			simple <- TRUE
			message("resorting to simple = ", simple)
		}
	}
	name.lab <- x@statistic.name
	z.lab <- "z-transformed p-value"
	graph <- function(X, Y, Xlab, Ylab){plot(x = X, y = Y, xlab = Xlab, ylab = Ylab, ...)}
	if(is.null(simple))
		graph(X = x@statistic, Y = x@zz, Xlab = name.lab, Ylab = z.lab)
	else if(simple)
		graph(Y = x@statistic, X = x@zz, Ylab = name.lab, Xlab = z.lab)
})

setClass("cvalues", representation("list", b = "numeric", b.name = "character"))
setValidity("cvalues", function(object)
{
	len.ok <- length(object) == length(object@b) && length(object@b.name) == 1
	cla.ok <- all(sapply(object, is, class2 = "cvalue"))
	ok <- len.ok && cla.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("[", signature(x = "cvalues", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@.Data <- as(x, "list")[i]
  x@b <- x@b[i]
  stopifnot(validObject(x))
  x
})
setMethod("plot", signature(x = "cvalues", y = "cvalues"), function(x, y, xlegend, ylegend, ylab = "number of sign calls", i, separate.legend = FALSE, xlab, ...)
{
	if(missing(xlab))
		xlab <- x@b.name
	if(!missing(i))
	{
		x <- x[i]
		y <- y[i]
	}
	stopifnot(length(x) == length(y) && all(x@b == y@b))
	is.stat <- function(object){all(sapply(x, is, "statistic.cvalue"))}
	if(is.stat(x) && is.stat(y))
	{
		sta <- function(object){sapply(object, function(cval){sum(cval@statistic != 0, na.rm = TRUE)})}
		x.pch <- 2
		y.pch <- 6
		xstat <- sta(x)
		ystat <- sta(y)
		B <- ifelse(x@b >= 0, x@b, 2/3)
		if(separate.legend) Mfrow()
		plot(x = B, y = xstat, xlab = xlab, ylab = ylab, log = "xy", ylim = range(xstat, ystat, na.rm = TRUE), pch = x.pch, ...)
		points(x = B, y = ystat, pch = y.pch)
		if(separate.legend) blank.plot()
		legend(legend = c(xlegend, ylegend), x = if(separate.legend) "left" else "topright", pch = c(x.pch, y.pch))
	}
	else
		stop("cannot plot it")
})

setClass("conditioningMerit", representation(nonancillarity = "Numeric", relevance = "Scalar", nsilence = "Numeric"))
setValidity("conditioningMerit", function(object)
{
	ok <- len.ok <- length(object@nonancillarity) == length(object@nsilence)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "conditioningMerit", y = "missing"), function(x, y, text.x, text.y, ...)
{
	plot(x = x@nsilence, y = x@nonancillarity, xlab = "number of affected features", ylab = "nonancillarity", log = "", ...)
	abline(h = x@relevance, col = "gray")
	if(missing(text.x)) text.x <- max(x@nsilence) / 7
	if(missing(text.y)) text.y <- x@relevance
	text(x = text.x, y = text.y, labels = "relevance")
})
setMethod("plot", signature(x = "conditioningMerit", y = "conditioningMerit"), function(x, y, text.x, text.y, ...)
{
	plot(x = x@nonancillarity, y = y@nonancillarity)
	abline(h = y@relevance, v = x@relevance, col = "gray")
	abline(a = 0, b = 1, col = "gray")
})

setClass("coverage", representation("numeric", n0 = "Scalar", nominal.rate = "Scalar", p0 = "Scalar", mean = "scalar", sd = "Scalar", mean0 = "scalar", sd0 = "Scalar", nulltype = "Scalar"))
setValidity("coverage", function(object)
{
	ok <- all(object <= object@n0)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "coverage", y = "missing"), function(x, y, ...)
{
	rate <- as(x, "numeric") / x@n0
	nulltype <- x@nulltype
	nulltype.name <- if(nulltype == 0) # theor
		"assumed null"
	else if(nulltype == 1) # MLE
		"estimated null"
	else if(nulltype == 2) # central matching
		"centrally matched null"
	else
		stop("nulltype not recognized")
	hist(rate, xlab = "actual coverage", ylab = "number of simulation trials", main = paste(100 * x@p0, "% null; ", nulltype.name, sep = ""), sub = paste("(", 100 * round(mean(rate) - x@nominal.rate, 3), "% coverage bias)", sep = ""))
})

setClass("riskPair", representation(risk1 = "numeric", risk0 = "numeric", FUN = "function", loss.FUN = "function")) # "Numeric" changed to "numeric" 090718
setValidity("riskPair", function(object)
{
	ok <- length(object@risk1) == length(object@risk0)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("plot", signature(x = "riskPair", y = "missing"), function(x, y, prob.name, nulltype.name, call.par = FALSE, mse, n0, p0, breaks = 50, ...) # called by "plot", signature(x = "probabilisticRisk", y = "missing")
{
	if(breaks >= length(x) + 1)
		breaks <- "Sturges"
	if(missing(mse))
		mse <- default(!is(x@loss.FUN, "cvalueError.FUN"), "mse")
	if(call.par)
		par(mfrow = c(2, 2))
	trans <- function(vec){if(mse) sqrt(vec) else vec}
	risk.name <- if(mse)
		"RMSE"
	else if(is(x@loss.FUN, "cvalueError.FUN"))
		x@loss.FUN@ann
	else
		"risk"
	digits <- 3
	sig <- function(x){paste(100 * signif(x, digits), "%", sep = "")}
	subplot <- function(X, Main, portion)
	{
		sub <- if(mse)
		{
			average <- stat(trans(X), FUN = median)
			weighted <- stat(X, FUN = mean) * portion
			av.text <- paste("med. RMSE: ", signif(average, digits), sep = "")
			weighted.text <- paste("wt. MSE: ", signif(weighted, digits), sep = "")
			paste(av.text, weighted.text, sep = "; ")
		}
		else
		{
			Sta <- Stats(trans(X))
			Q1 <- Sta["Q1"]
			average <- Q2 <- Sta["Q2"]
			Q3 <- Sta["Q3"]
			mean.abs <- Sta["mean.abs"]
			paste("abs.: ", sig(mean.abs), " (", sig(Q1), ", ", sig(Q2), ", ", sig(Q3), ")", sep = "")
#			paste("abs.: ", sig(mean.abs), " (Q1: ", sig(Q1), ", Q2: ", sig(Q2), ", Q3: ", sig(Q3), ")", sep = "")
		}
		hist(x = trans(X), xlab = paste(risk.name, " of ", prob.name, sep = ""), ylab = "number of trials", main = paste(Main, " features", sep = ""), sub = sub, breaks = breaks, ...)
		abline(v = average, col = "gray")
		legend(legend = paste("(", nulltype.name, ")", sep = ""), bty = "n", x = if(mse) "topright" else "top")
		sig(mean.abs) # data.frame(status = if(affected) "affected" else "unaffected")
	}
	get.Main <- function(pMain, tex)
	{
		paste(100 * signif(pMain, digits), "% ", tex, sep = "")
	}
	affected <- subplot(X = x@risk1, Main = get.Main(1 - p0, "affected"), portion = 1 - p0)
	unaffected <- subplot(X = x@risk0, Main = get.Main(p0, "unaffected"), portion = p0)
#	invisible(rbind(datf1, datf2))
	invisible(data.frame(p1 = sig(1 - as(p0, "numeric")), affected = affected, unaffected = unaffected))
})
setMethod("length", signature(x = "riskPair"), function(x)
{
	length(x@risk1)
})

setClass("Scalar.factor", representation("Scalar")) # used in ran.cvalue
setClass("random.scalar", representation(fixed = "scalar", random = "function"))
setClassUnion("extended.scalar", c("scalar", "random.scalar"))

setClass("cvalueError.FUN", representation("function", ann = "character", min_congruity_accept = "Scalar", max_congruity_accept = "Scalar"))
setValidity("cvalueError.FUN", function(object)
{
	lt1 <- function(x){length(x) == 1 && !is.na(x) && x <= 1.0001}
	range.ok <- lt1(object@min_congruity_accept) && lt1(object@max_congruity_accept) && object@min_congruity_accept <= object@max_congruity_accept
	ann.ok <- length(object@ann) == 1 && !is.na(object@ann)
	ok <- range.ok && ann.ok
	if(!ok)
	{
		printInvalid(object)
		browser()
	}
	ok
})


setClass("probabilisticRisk", representation(cvalue.risk = "riskPair", fdr.risk = "riskPair", n0 = "Scalar", p0 = "Scalar", mean = "extended.scalar", sd = "extended.scalar", mean0 = "extended.scalar", sd0 = "extended.scalar", nulltype = "Scalar", loss.FUN = "function"))
setValidity("probabilisticRisk", function(object)
{
	p0.ok <- object@p0 <= 1
	len.ok <- try(length(object@cvalue.risk) == length(object@fdr.risk))
	ok <- !is(len.ok, "try-error") && (len.ok && p0.ok)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("length", signature(x = "probabilisticRisk"), function(x)
{
	length(x@cvalue.risk)
})
setMethod("plot", signature(x = "probabilisticRisk", y = "missing"), function(x, y, call.par = TRUE, plot.fdr.risk = default(FALSE, "plot.fdr.risk"), ...)
{
	nulltype <- x@nulltype
	nulltype.name <- if(nulltype == 0) # theor
		"assumed null"
	else if(nulltype == 1) # MLE
		"estimated null"
	else if(nulltype == 2) # central matching
		"centrally matched null"
	else
		stop("nulltype not recognized")
	subplot <- function(X, prob.name)
	{
		plot(x = X, prob.name = prob.name, nulltype.name = nulltype.name, mse = identical(x@loss.FUN, squaredError), n0 = x@n0, p0 = x@p0, ...)
	}
	if(call.par)
		par(mfrow = c(2, 2))
	datf <- subplot(X = x@cvalue.risk, prob.name = "certainty level") # 	invisible(data.frame(p1 = 1 - p0, affected = affected, unaffected = unaffected))
	datf$nulltype <- nulltype
	if(plot.fdr.risk)
		subplot(X = x@fdr.risk, prob.name = "local FDR")
	invisible(datf)
})
setAs(from = "probabilisticRisk", to = "data.frame", function(from)
{
#	data.frame(p1 = 1 - from@p0, )
})

setClass("dataGenerator", representation("function", CDF0 = "CDF", pvalue.fun = "function", zz = "numeric"))
setMethod("plot", signature(x = "dataGenerator", y = "missing"), function(x, y, main = "dataGenerator", congruity.range.accept = numeric(0), normal = logical(0), call.par = length(normal) == 0, ...)
{
#	plot(x@CDF0, ...)
	zcdf0 <- function(zz) # plot(rfun0, main = paste(n, "features"), xlim = c(-3, 3))
	{
		qnorm(x@CDF0(zz))
	}
	zz <- seq(-3, 3, length.out = 1000)
	if(call.par)
		Mfrow()
	graph <- function(norm)
	{
		xy <- function(xy){if(norm) xy else pnorm(xy)}
		plot(xy(zz), xy(sapply(zz, zcdf0)), main = main, ...)
		abline(a = 0, b = 1, col = "gray")
		if(length(congruity.range.accept) == 2)
		{
			vh <- if(norm)
				qnorm(congruity.range.accept)
			else
				congruity.range.accept
			abline(v = vh, h = vh, col = "gray")
		}
	}
	if(length(normal) == 0)
	{
		graph(FALSE)
		graph(TRUE)
	}
	else
		graph(normal)
})

setClassUnion("marginal.cvalue", c("unconditional.cvalue", "NULL"))
setClass("probabilisticRisks", representation("list", file = "character", marginal.cval = "marginal.cvalue"))
setValidity("probabilisticRisks", function(object)
{
	cla.ok <- all(sapply(object, is, class2 = "probabilisticRisk"))
	len.ok <- length(object) > 0 && length(object@file) %in% c(0, 1)
	ok <- cla.ok && len.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("stripchart", signature(x = "probabilisticRisks"), function(x, xlim, cex = 2, call.legend = default(TRUE, "call.legend"), call.par = call.legend, main = class(x), sub = "", all.Q = default(FALSE, "all.Q"), ...)
{
	assert.is(x[[1]]@loss.FUN, "cvalueError.FUN")
	if(call.par)
		par(mfrow = c(2, 2))
	min.accept <- min_congruity_accept(x)
	max.accept <- max_congruity_accept(x)
	ann <- annotation(x)
	nam <- if(is.null(names(x)))
		make.names(1:length(x))
	else
		names(x)
	if(missing(xlim))
	{
		xlim <- max(abs(range(as.numeric(sapply(x, function(risk){c(risk@cvalue.risk@risk1, risk@cvalue.risk@risk0)})), na.rm = TRUE)))
		xlim <- c(-xlim, xlim)
	}
	graph <- function(add, Stats.name, risk.name, col, pch)
	{
		summarize <- function(risk)
		{
			crisk <- risk@cvalue.risk
			vec <- Stats(slot(crisk, name = risk.name))[Stats.name]
			if(!all(is.finite(vec)))
			{ message("bad summarize"); browser()}
			vec
		}
		lis <- lapply(x, summarize)
		names(lis) <- nam
		stripchart(x = lis, add = add, col = col, pch = pch, xlim = xlim, xlab = paste(ann, "of certainty level"), cex = cex, ...)
		title(main = main, sub = sub)
	}
	pch0 <- 1
	pch1 <- 4
	Q.col <- "gray"
	mean.abs.col <- "black"
	Q.name <- if(all.Q) c("Q1", "Q2", "Q3") else "Q2"
	mean.abs.name <- "mean.abs"
	graph(add = FALSE, Stats.name = Q.name, risk.name = "risk0", col = Q.col, pch = pch0)
	graph(add = TRUE, Stats.name = Q.name, risk.name = "risk1", col = Q.col, pch = pch1)
	graph(add = TRUE, Stats.name = mean.abs.name, risk.name = "risk0", col = mean.abs.col, pch = pch0)
	graph(add = TRUE, Stats.name = mean.abs.name, risk.name = "risk1", col = mean.abs.col, pch = pch1)
	if(call.legend)
	{
		blank.plot()
		Q.leg <- if(all.Q) "quartiles;" else "median;"
		mean.abs.leg <- "mean abs.;"
		leg0 <- "zero mean"
		leg1 <- "nonzero mean"
		legend(x = "left", legend = paste(c(Q.leg, Q.leg, mean.abs.leg, mean.abs.leg), c(leg0, leg1, leg0, leg1)), col = c(Q.col, Q.col, mean.abs.col, mean.abs.col), pch = c(pch0, pch1, pch0, pch1))
	}
})
setMethod("plot", signature(x = "probabilisticRisks", y = "missing"), function(x, y, file, plot.fdr.risk = default(FALSE, "plot.fdr.risk"), ...)
{
	if(missing(file))
		stripchart(x = x, ...)
	else
	{
		message("plotting to ", file, " on ", date())
		Pdf(file = paste(file, sep = ".")) # width = big.width, height = big.height, pointsize = big.pointsize
		lis <- lapply(1:length(x), function(i)
		{
			call.par <- if(plot.fdr.risk)
				TRUE
			else
				i %% 2 == 1
			plot(x[[i]], call.par = call.par, plot.fdr.risk = plot.fdr.risk, ...)
		})
		dev.off()
		datf0 <- do.call("rbind", lis) # p1, affected, unaffected, nulltype
		stopifnot(nrow(datf0) == length(x))
		npair <- nrow(datf0) / 2
		stopifnot(npair == floor(npair))
		lis1 <- lapply(1:npair, function(i)
		{
			odd <- datf0[i * 2 - 1, ]
			even <- datf0[i * 2, ]
			stopifnot(odd$p1 == even$p1)
			assumed <- if(odd$nulltype == 0)
				odd
			else if(even$nulltype == 0)
				even
			else
				NULL
			estimated <- if(odd$nulltype > 0)
				odd
			else if(even$nulltype > 0)
				even
			else
				NULL
			if(is.data.frame(assumed) && is.data.frame(estimated))
			{
				dat <- try(data.frame(p1 = odd$p1, affected.estimated = estimated$affected, affected.assumed = assumed$affected, unaffected.estimated = estimated$unaffected, unaffected.assumed = assumed$unaffected))
				if(is(dat, "try-error"))
				{ message("dat err"); print(sapply(list(p1 = odd$p1, affected.estimated = estimated$affected, affected.assumed = assumed$affected, unaffected.estimated = estimated$unaffected, unaffected.assumed = assumed$unaffected), length)); browser()}
				dat
			}
			else
			{  message("unexpected structure; cannot create data frame"); browser()}
		})
		datf1 <- do.call("rbind", lis1)
		invisible(datf1)
	}
})
setMethod("[", signature(x = "probabilisticRisks", i = "ANY", j = "missing"), function(x, i, j, drop)
{
  x@.Data <- as(x, "list")[i]
  stopifnot(validObject(x))
  x
})



setMethod("print", signature(x = "empiricalNull"), function(x)
{
	rou <- function(...){round(..., digits = 3)}
	datf <- data.frame(Mean = rou(Mean(x, nulltype = 1)), Sd = rou(Sd(x, nulltype = 1)), relevance = rou(halfAlphaI(x)))
	print(datf)
})
setMethod("plot", signature(x = "numeric", y = "empiricalNull"), function(x, y, main = y@CDF0@param.name, cval, acval, ...)
{
	par(mfrow = c(2, 2))
#	plot(y@CDF0)
	Main <- main
	plot.fun <- function(fun, ylab, main)
	{
		if(missing(main))
			main <- Main
		plot(x = x, y = fun, ylab = ylab, main = main, xlab = "z-transform of certainty level", ...)
		abline(v = 0, col = "gray")
	}
	cval <- pnorm(x)
	names(cval) <- names(x)
	Y <- sapply(x, as(y@CDF0, "function"))
	plot(x = cval, y = Y, xlab = "nominal certainty level", ylab = "adjusted certainty level", main = "certainty level adjustments", ylim = c(0, 1), type = "l")
	abline(0, 1, col = "gray")
	abline(v = 0.5, h = 0.5, col = "gray")
	if(!missing(cval) && !missing(acval))
	{
		assert.is(cval, "cvalue")
		assert.is(acval, "adjusted.cvalue")
		plot(x = cval, y = acval)
	}
	else if(missing(cval) && missing(acval)) # original implementation
	{
		plot.fun(fun = y@PDF, ylab = "density estimate")
		plot.fun(fun = y@PDF0, ylab = "null density estimate")
		hist(lfdr(y), xlab = "local FDR estimate", main = paste(round(y@p0, 3) * 100, "% null", sep = ""))
	#	plot.fun(fun = as(y, "function"), ylab = "local FDR estimate", main = paste(round(y@p0, 3) * 100, "% null", sep = "")) # this conflicts with lfdr for an unknown reason
	}
	else if(!missing(cval) && missing(acval)) # added 100126
	{
		message(class(cval), " ", class(x), " ", class(y))
		if(is(cval, "cvalue"))
			plot(x = cval, y = y, call.par = FALSE, ...) # recurse = FALSE,
		else if(is.numeric(cval))
		{
			stop('consider "plot", signature(x = "cvalue", y = "empiricalNull") with simple = TRUE')
			fdr <- lfdr(y)
			if(length(cval) != length(fdr))
			{
				nam <- intersect(names(cval), names(fdr))
				nam.ok <- (is.character(nam) && length(nam) > 0)
				if(!nam.ok)
				{ message("bad nam in plot method"); browser()}
				cval <- cval[nam]
				fdr <- fdr[nam]
			}
			graph <- function(p, plab)
			{
				plot(x = p, y = fdr, xlab = plab, ylab = "local false discovery rate", ...)
			}
			graph(p = x, plab = "z")
			graph(p = cval, plab = "p value")
			message('consider "plot", signature(x = "cvalue", y = "empiricalNull") with simple = TRUE')
		}
		else
			stop("plot method not implemented for arguments given!!!!!!!!")
	}
	else
		stop("plot method not implemented for arguments given...")
})
setMethod("plot", signature(x = "empiricalNull", y = "missing"), function(x, y, ...)
{
	xx <- seq(x@CDF0@min_param, x@CDF0@max_param, length.out = 100)
	if(!all(is.finite(xx)))
	{ message("try plot(x = [numeric], y = x, ...)"); browser()}
	plot(x = xx, y = x, ...)
})
setMethod("plot", signature(x = "empiricalNull", y = "cvalue"), function(x, y, ...)
{
	plot(y, x, ...)
})
setMethod("plot", signature(x = "cvalue", y = "empiricalNull"), function(x, y, statistic.range, simple = FALSE, call.ecdf = default(!simple, "call.ecdf"), lwd = 3, call.par = TRUE, adjust = !simple, ...)
{
	if(call.par)
		par(mfrow = c(2, 2))
	nam <- intersect(names(x), names(y))
	plot(pnorm(x@zz)[nam], lfdr(y)[nam], xlab = "certainty level", ylab = "estimated local false discovery rate")
	if(missing(statistic.range)) # statistic.range is like mu.range of odd.s
	{
		if(call.ecdf)
		{
			ecdf.col <- "gray"
			null.col <- "black"
			cdf0 <- y@CDF0
			X <- pnorm(seq(cdf0@min_param, cdf0@max_param, length.out = 1e3))
			plot(x = X, y = sapply(X, function(p){cdf0(qnorm(p))}), xlab = "nominal certainty level", ylab = "cumulative probability estimate", main = "estimated certainty level CDFs", col = null.col, type = "l", lwd = lwd, ...)
			abline(v = 0.5, h = 0.5, col = "gray")
			abline(a = 0, b = 1, col = "gray")
			xcdf <- function(q){ecdf(congruity(x))(q)} # ECDF(congruity(x))
			assert.is(xcdf, "function")
			plot(x = X, y = xcdf, col = ecdf.col, add = TRUE, lwd = lwd, ...)
#			plot(x = X, y = xcdf, col = ecdf.col, ...)
			legend(x = "topleft", legend = c("null CDF estimate", "empirical CDF"), col = c(null.col, ecdf.col), bty = "n", lty = 1, lwd = lwd)
		}
		else
		{
			plot.pair <- function(cval, main)
			{
				subplot <- function(X, xlab)
				{
					Plot(x = X, y = lfdr(y), xlab = xlab, ylab = "local FDR estimate", main = main, ...)
				}
				subplot(X = cval@zz, xlab = "z-transform of certainty level")
				subplot(X = congruity(cval), xlab = "congruity")
			}
			plot.pair(cval = x, main = "nominal certainty level")
			if(adjust)
				plot.pair(cval = adjusted.cvalue(x, plot = 0), main = "adjusted certainty level")
		}
	}
	else
	{
		stop("not yet implemented for statistic.range")
	}
})
setMethod("Plot", signature(x = "cvalue", y = "empiricalNull"), function(x, y, ...)
{
	get.con <- function(...){confidence(object = x, ...)}
	alt.con <- get.con(less.than = TRUE)
	graph <- function(less.than, fun, ...)
	{
		mixture.con <- get.con(mass0 = y, less.than = less.than)
		mass0 <- Mass0(y)#XXX|:no visible global function definition for ÔMass0Õ: in interval.s
		nam <- intersect(names(alt.con), names(mixture.con))
		nam <- intersect(nam, names(mass0))
		alt.con <- alt.con[nam]
		mixture.con <- mixture.con[nam]
		mass0 <- mass0[nam]
		ord <- order(alt.con)
		fun(x = alt.con[ord], y = mixture.con[ord], ...)
	}
	less.lty <- 1
	greater.lty <- 2
	less.col <- "black"
	greater.col <- "gray"
	graph(less.than = TRUE, fun = plot, type = "l", lty = less.lty, xlab = "conditional confidence that mean < 0", ylab = "marginal confidence...", col = less.col, ...)
	graph(less.than = FALSE, fun = lines, lty = greater.lty, col = greater.col)
	legend(x = "top", legend = c("... that mean < 0", "... that mean > 0"), lty = c(less.lty, greater.lty), col = c(less.col, greater.col))
})
setMethod("plot", signature(x = "empiricalNull", y = "empiricalNull"), function(x, y, x.i = 1:length(x.fdr), y.i = 1:length(y.fdr), ...)
{
	x.fdr <- lfdr(x)
	y.fdr <- lfdr(y)
	X <- x.fdr[x.i]
	Y <- y.fdr[y.i]
	if(length(X) != length(Y))
	{ message("bad length!"); browser()}
	plot(X, Y, main = "local false discovery rates", ...)
})

setMethod("names", signature(x = "empiricalNull"), function(x)
{
	names(lfdr(x))
})
setMethod("length", signature(x = "empiricalNull"), function(x)
{
	length(lfdr(x))
})
setAs(from = "empiricalNull", to = "function", function(from)
{
	warning("this conflicts with lfdr for an unknown reason")
	function(q)
	{
		fdr <- from@p0 * from@PDF0(q) / from@PDF(q)
		ifelse(fdr <= 1, fdr, 1)
	}
})
setAs(from = "empiricalNull", to = "CDF", function(from)
{
#	new.CDF(object = object, min_param = from@CDF0@min_param, max_param = from@CDF0@max_param, param.name = from@CDF0@param.name, type = from@CDF0@type)
	from@CDF0
})


# near end of file:

#source(file = "congruity.s") # functions











