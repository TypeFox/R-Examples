# MDL.s created by David Bickel on 14 July 2010.

E <- function(x, P, check.P = FALSE, ...)
{
	assert.are(list(x, P), "numeric")
	if(length(x) == 1)
		x <- rep(x, length(P))
	else if(length(P) == 1)
		P <- rep(P, length(x))
	tol <- 1e-3
	ok <- if(check.P)
		abs(sum(P, na.rm = TRUE) - 1) <= tol
	else
		is.prob(P = P, tolerance = tol, ...)
	stopifnot(ok)
	scalar(sum(x * P, ...))
}
entropy <- function(P, base = exp(1)) # see divergence
{
	stop("This should be patterned after divergence.")
	vec <- ifelse(P == 0, 0, E(x = -logb(P, base = base), P = P))
	ok <- length(vec) == length(P) && all(vec >= 0)
	stopifnot(ok)
	vec
}
Entropy <- function(P, base = exp(1), ...)
{
	base ^ entropy(P = P, base = base, ...)#XXX|check.P ???:unused argument (check.P = FALSE), I removed it:check.P = FALSE,  
}
prob.density <- function(object, x, base = exp(1), scalar.density = TRUE, ...)
{
	assert.is(x, "numeric")
	assert.is(base, "numeric")
	stopifnot(length(base) == 1)
	dens <- if(is(object, "univariateDistribution"))
	{
		dfun <- d(object)
		assert.is(dfun, "function")
		prob.density(object = dfun, x = x, base = base, scalar.density = scalar.density, ...)
	}
	else if(is(object, "DensityFunNML"))
	{
		dfun <- function(X, base = exp(1), log = TRUE)
		{
			Dens <- object(X)
			assert.is(Dens, "numeric")
			if(log)
				logb(Dens, base = base)
			else
				Dens
		}
		if(is(dfun, "DensityFunNML"))
			stop("recursive DensityFunNML")
		prob.density(object = dfun, x = x, base = base, scalar.density = scalar.density, ...)
	}
	else if(is(object, "function"))
	{
		densities <- object(x, log = TRUE, ...)
		ln.dens <- if(scalar.density)
			sum(densities) # IID
		else
			densities
		if(base == exp(1))
			ln.dens
		else
			new.log(ln.dens, base = base)
	}
	else
		stop("bad prob.density")
	if(scalar.density)
		scalar(dens)
	else
	{
		names(dens) <- names(x)
		dens
	}
}

Df <- function(x)
{
	if(is(x, "xprnSet"))
		ncol(exprs(x)) - 1
	else if(is(x, "xprnSetPair"))
		ncol(exprs(x@x)) + ncol(exprs(x@y)) - 2
	else if(is(x, "matrix"))
		ncol(x)-1 #added by marta
	else
		df(x)
}

new.Family <- function(Distr.fun, unknown.param.name, known.param, discrete, lower.unknown.param, upper.unknown.param)
{
	assert.is(unknown.param.name, "character")
	new("Family", Distr.fun = Distr.fun, unknown.param.name = unknown.param.name, known.param = known.param, discrete = discrete, lower.unknown.param = scalar(lower.unknown.param), upper.unknown.param = scalar(upper.unknown.param))
}

setMethod("is.nothing", signature(object = "Family"), function(object)
{
	nam <- slotNames(object)
	nam <- nam[nam != "Distr.fun"]
	are.nothing(lapply(nam, slot, object = object))
})
blank.Family <- function()
{
	object <- new.Family(Distr.fun = identity, unknown.param.name = character(0), known.param = list(), discrete = logical(0), lower.unknown.param = scalar(numeric(0)), upper.unknown.param = scalar(numeric(0)))
	stopifnot(is.nothing(object))
	object
}
new.Families <- function(object){new("Families", object)}
####====================================
Family <- function(object, unknown.param.name, known.param, discrete, lower.unknown.param, upper.unknown.param, lower = 0, ...)
{
	arglis <- list(...)
	if(is(object, "function") && length(arglis) == 0)
	{
		if(is(object, "likFun"))
			stop("object does not return univariateDistribution object")
		assert.is(unknown.param.name, "character")
		new.Family(Distr.fun = object, unknown.param.name = unknown.param.name, known.param = known.param, discrete = discrete, lower.unknown.param = lower.unknown.param, upper.unknown.param = upper.unknown.param)
	}
	else if(is(object, "xprnSetObject") && missing(unknown.param.name) && missing(known.param) && missing(discrete) && missing(lower.unknown.param) && missing(upper.unknown.param))
	{
		extendedFamily.absTd(df = Df(object), lower = lower, ...)
	}
	else
	{	message("cannot generate Family from arguments given"); browser()}
}
new.extendedFamily <- extendedFamily <- function(object, mle.fun)
{
	assert.is(object, "Family")
	new("extendedFamily", object, mle.fun = mle.fun)
}
unknownParam <- function(object)
{
	unknown.param <- object
	if(is(unknown.param, "resampled.posteriorP0"))
		unknown.param <- unknown.param@original
	if(is(unknown.param, "posteriorP0"))
		unknown.param <- unknown.param@estimate
	if(is(unknown.param, "nmle") || is(unknown.param, "posteriorP0.error"))#XXX|:changing from mle to nmle
		unknown.param <- unknown.param@unknown.param
	if(is(unknown.param, "numeric"))
		unknown.param <- as.list(unknown.param) # good when unknown.param was returned by an mle.fun slot?
	assert.is(unknown.param, "list") # added 110406
	if(is(unknown.param, "unknownParams"))
		stop("bad unknown.param")
	unknown.param
}
new.unknownParams <- function(object)
{
	if(is(object, "unknownParams"))
		object
	else
	{
		lis <- lapply(object, unknownParam)
		new("unknownParams", lis)
	}
}
param.list <- function(object, unknown.param)
{
	assert.is(object, "Family")
	unknown.param.name <- object@unknown.param.name
	known.param <- object@known.param
	unknown.param <- unknownParam(object = unknown.param)
	if(is.null(names(unknown.param)))
	{
		if(length(unknown.param) == 1 && length(unknown.param.name) == 1)
			names(unknown.param) <- unknown.param.name
		else
			stop("ambiguous parameter arguments in likelihood function")
	}
	blank <- as.list(unknown.param.name)
	names(blank) <- unknown.param.name
	stopifnot(all(unknown.param.name %in% names(unknown.param)))
#	unknown.param <- unknown.param[unknown.param.name] # good when unknown.param was returned by an mle.fun slot?
	if(!sameNames(blank, unknown.param, order.sensitive = FALSE))
	{ message("unknown.param does not match unknown.param.name"); browser()}
	assert.are(list(known.param, unknown.param), "list")
	c(known.param, unknown.param)
}
Distr <- function(object, unknown.param)
{
	assert.is(object, "Family")
	param <- param.list(object = object, unknown.param = unknown.param)
	dis <- do.call(object@Distr.fun, param)
	assert.is(dis, "univariateDistribution")
	if(object@discrete)
		assert.is(dis, "discreteDistribution")
	else
		assert.is(dis, "abscontDistribution")
	dis
}

statistic.Td <- function(x, y, cv, size, fam = NULL, na.rm = FALSE, verbose = FALSE)
{
	if(na.rm && verbose)
		warning("statistic.Td is not implemented for missing data; 'na.rm = TRUE' will be ignored")
	if(!missing(fam))
		stopifnot(identical(fam@Distr.fun, Td))
	df.ok <- function(df)
	{
		is.null(fam) || fam@known.param$df %in% df #stopifnot(size.ok(size = sum(is.finite(x))))
	}
	get.num <- function(num){num[is.finite(num)]}
	stat <- if(missing(y))
	{
		if(missing(x))
		{
			stopifnot(df.ok(df = size - 1))
			new.statistic(sqrt(size) / cv)
		}
		else if(missing(cv) && missing(size))
		{
			if(!df.ok(df = sum(is.finite(x)) - 1))
			{ message("degree of freedom problem 1"); browser()}
			x <- get.num(num = x)
			statistic.Td(cv = sd(x) / mean(x), size = length(x), na.rm = na.rm)
		}
	}
	else if(is(x, "numeric") && is(y, "numeric"))
	{
		if(!df.ok(df = sum(is.finite(c(x, y))) - 2))
		{ message("degree of freedom problem 2"); browser()}
		x <- get.num(num = x)
		y <- get.num(num = y)
		new.statistic(t_statistic(x = x, y = y)) # new.statistic(t.test(x = x, y = y, var.equal = TRUE)$stat) is too slow
	}
	else
		stop("cannot compute t statistic from the arguments given")
	stopifnot(length(stat) == 1)
	stat
}
statistic.absTd <- function(..., fam)
{
	tstat <- if(missing(fam))
		statistic.Td(...)
	else
	{
		if(!sameFunction(fam@Distr.fun, absTd))
		{ message("bad fam"); browser()}
		df <- fam@known.param$df
		statistic.Td(..., fam = Family.Td(df = df))
	}
	abs(tstat)
}
new.statistic <- function(object)
{
	if(is(object, "statistic"))
		object
	else
		new("statistic", object)
}
#XXX|:change name: from statistic to nstatistic
nstatistic <- function(object = NULL, x, y = numeric(0), x.converges = logical(0), size = numeric(0), na.rm = FALSE)
{
	if(missing(x))
	{
		if(is(object, "xprnSetObject"))
		{
			x <- object
			object <- NULL
		}
		else
			stop("x is missing, and statistic could not determine it from object")
	}
	if(is.null(object))
	{
		if(is(x, "xprnSetObject"))
			object <- extendedFamily.absTd(df = Df(x))
		else
			stop("object is missing, and statistic could not determine the Family from x")
	}
	assert.is(object, c("Family", "Families"))
	reduced <- if(is(x, "xprnSetObject"))
	{
		get.stat <- function(x1, y1, fam)
		{
			assert.is(fam, "Family")
			stat(x = x1, y = y1, FUN = nstatistic, object = fam, x.converges = x.converges, size = size, na.rm = na.rm) # should this call statistic instead of stat?
		}
		if(is(object, "Family"))
			get.stat(x1 = x, y1 = y, fam = object)
		else if(is(object, "Families") && ncomparison(x) == length(object))
		{
			reduced.x <- sapply(1:length(object), function(i)
			{
				x1 <- x[i, ]
				y1 <- if(is.nothing(y))
					y
				else if(is(y, "xprnSetObject"))
					y[i, ]
				else
					stop("cannot reduce the expression data")
				get.stat(x1 = x1, y1 = y1, fam = object[[i]])
			})
			names(reduced.x) <- rownames(x)
			reduced.x
		}
		else
			stop("statistic problem")
	}
	else if(is(x, "numeric") && is(object, "Family"))
	{
		if(is.nothing(x.converges) && is.nothing(size) && length(x) >= 2) # 100802
		{
			if(is(x, "statistic"))
				stop("statistic of a statistic?")
			get_statistic <- function(fun)
			{
				if(is.nothing(y))
					fun(x = x, fam = object, na.rm = na.rm)
				else
					fun(x = x, y = y, fam = object, na.rm = na.rm)
			}
			fun <- if(is.Family.Td(object))
				statistic.Td
			else if(is.Family.absTd(object))
				statistic.absTd
			else
			{ message("cannot reduce the data to a statistic, perhaps because not yet implemented for the family"); browser()}
			sta <- get_statistic(fun = fun)
			if(!is.finite(sta))
			{ message("non-finite statistic"); browser()}
			sta
		}
		else if(!is.nothing(x.converges) && is.nothing(y)) # before 100802
		{
			if(is(x, "statistic") && x.converges)
				stop("double transform statistic")
			y <- if(x.converges && is.Family.Td(object) && !is.nothing(size))
			{
				if(is.nothing(size))
					size <- object@known.param$df + 1
				statistic.Td(cv = x, size = size)
			}
			else if(!x.converges) #  && is.nothing(size)
				x
			else
				stop("conversion statistic not implemented for that object@Distr.fun")
			assert.is(y, "numeric")
			y
		}
		else
			stop("statistic error")
	} # end if
	else
		stop("cannot compute statistic for this class")
	new.statistic(object = reduced)
} # end statistic
#
#####====================================
new.finiteFamily <- function(object, unknown.param.set)
{
	object@lower.unknown.param <- object@upper.unknown.param <- scalar(numeric(0))
	new("finiteFamily", object, unknown.param.set = unknown.param.set)
}
finiteFamily <- function(..., unknown.param.set)
{
	new.finiteFamily(Family(...), unknown.param.set = unknown.param.set)
}

new.likFun <- function(object, Distr.fun, unknown.param.name, known.param, x, base, weight = blank.Weight(), arglis)
{
	assert.is(unknown.param.name, "character")
	new("likFun", object, Distr.fun = Distr.fun, unknown.param.name = unknown.param.name, known.param = known.param, x = x, base = Scalar(base), weight = weight, arglis = arglis)
}
likFun <- function(object, unknown.param.name = character(0), x = numeric(0), known.param = list(), base = exp(1), discrete = logical(0), lower.unknown.param = numeric(0), upper.unknown.param = numeric(0), param0 = list(), ...) # The returned "arglis" slot is up to date with the argument list of likFun as of 7 September 2010.
{
	arglis <- list(...)
	if(is(object, "DensityFunNML") && is.nothing(x) && is.nothing(param0))
	{
		warning("not really a pseudo-likelihood function, but this is used for plotting purposes")
		lf <- likFun(object@family, x = as.numeric(NA), base = base)
		lf@.Data <- function(x, x.converges = FALSE)
		{
			Lik <- do.call(object, c(x = x, x.converges = x.converges, arglis, ...))
			logb(Lik, base = base)
		}
		lf
	}
	else if(is(object, "Family") && is.nothing(unknown.param.name) && is.nothing(known.param) && is.nothing(discrete) && is.nothing(lower.unknown.param) && is.nothing(upper.unknown.param) && length(arglis) == 0 && is(x, "numeric"))
	{
		fun <- function(...)
		{
			if(!all(names(param) %in% argnames(object@Distr.fun)))
			{ message("known.param or unknown.param does not match args(object@Distr.fun)"); browser()}
			unknown.param <- list(...)
			if(length(unknown.param) == 1 && is.null(names(unknown.param)))
			{
				if(length(unknown.param) == length(object@unknown.param.name))
					names(unknown.param) <- object@unknown.param.name
				else
				{ message("bad object@unknown.param.name?"); browser()}
			}
			blank <- as.list(object@unknown.param.name)
			names(blank) <- object@unknown.param.name
			if(!sameNames(blank, unknown.param, order.sensitive = FALSE))
			{ stop("unknown.param does not match object@unknown.param.name"); message("unknown.param does not match object@unknown.param.name"); print(unknown.param); browser()}
			if(length(unknown.param) == 1 && !((length(object@lower.unknown.param) == 0 || unknown.param >= object@lower.unknown.param) && (length(object@upper.unknown.param) == 0 || unknown.param <= object@upper.unknown.param)))
			{ message("likelihood function evaluated outside of domain"); browser()}
			dis <- Distr(object = object, unknown.param = unknown.param)
			assert.is(dis, "univariateDistribution")
			prob.density(object = dis, x = x, base = base)
		}
		lf <- if(is.nothing(param0))
			fun
		else if(is.list(param0))
			function(...){fun(...) - do.call(fun, param0)}
		else
			stop("bad param0")
		arglis <- c(list(object = object, unknown.param.name = unknown.param.name, x = x, known.param = known.param, base = base, discrete = discrete, lower.unknown.param = lower.unknown.param, upper.unknown.param = upper.unknown.param, param0 = param0), arglis) # enables recomputing likelihood with a changed argument, as required by "weighted.likFun"
		new.likFun(object = lf, Distr.fun = object@Distr.fun, unknown.param.name = object@unknown.param.name, known.param = object@known.param, x = x, base = base, arglis = arglis)
	}
	else if(is(object, "function") && length(arglis) == 0)
	{
		if(is(object, "likFun"))
			stop("object does not return univariateDistribution object")
		if(is.nothing(known.param))
			known.param <- list()
		assert.is(x, "numeric")
		assert.is(known.param, "list")
		assert.is(unknown.param.name, "character")
		likFun(object = Family(object = object, unknown.param.name = unknown.param.name, known.param = known.param, discrete = discrete, lower.unknown.param = lower.unknown.param, upper.unknown.param = upper.unknown.param), x = x, base = base, param0 = param0)
	}
	else
	{ message("bad likFun args"); browser()}
} # end likFun
new.Weight <- function(focus.weight = numeric(0), incidental.weight = numeric(0), lower.bound = numeric(0), upper.bound = numeric(0), heavy.focus = FALSE)
{
	if(!(is.nothing(focus.weight) && is.nothing(incidental.weight)))
	{
		tol <- 1e-3
		if(is.nothing(lower.bound))
			lower.bound <- 1 - tol
		if(is.nothing(upper.bound))
			upper.bound <- lower.bound + 2 * tol
		weight.ok <- function(focus.w, incidental.w){!heavy.focus || all(focus.w >= incidental.w)}
		if(!is.nothing(incidental.weight) && is.nothing(focus.weight))
		{
			assert.is(incidental.weight, "numeric")
			focus.weight <- 1
			if(!weight.ok(focus.w = focus.weight, incidental.w = incidental.weight))
			{ message("problem with numeric incidental.weight"); browser()}
			unnormalized.total <- sum(focus.weight, incidental.weight)
			normalized.total <- if(is.nothing(lower.bound) || is.nothing(upper.bound))
				1
			else if(!is.nothing(lower.bound) || !is.nothing(upper.bound))
				mean(c(lower.bound, upper.bound))
			else
				unnormalized.total
			norm <- normalized.total / unnormalized.total
			focus.weight <- norm * focus.weight
			incidental.weight <- norm * incidental.weight
		}
		if(!weight.ok(focus.w = focus.weight, incidental.w = incidental.weight))
		{ message("problem with weights"); browser()}
	}
	new("Weight", focus.weight = scalar(focus.weight), incidental.weight = incidental.weight, lower.bound = scalar(lower.bound), upper.bound = scalar(upper.bound))
}
setMethod("is.nothing", signature(object = "Weight"), function(object)
{
	is.nothing(object@focus.weight) && is.nothing(object@incidental.weight)
})
blank.Weight <- function()
{
	object <- new.Weight()
	stopifnot(is.nothing(object))
	object
}
new.Weights <- function(object){new("Weights", object)}
#===================
ncomparison <- function(object)
{
	if(is(object, "Weight"))
		1 + length(object@incidental.weight)
	else if(is(object, "Weights") || is(object, "statistic"))
		length(object)
	else if(is(object, "xprnSetObject"))
		length(names(object))
	else
	{	message("cannot compute ncomparison"); print(object); print(class(object)); browser()}
}

MWLE <- function(object, focus.x, incidental.x, focus.weight, incidental.weight, ...)
{# maximum weighted likelihood estimate

	assert.is(is(object, "Family"))
	assert.are(list(focus.x, incidental.x, focus.weight, incidental.weight), "numeric")
	stopifnot(length(focus.x) == length(focus.weight))
	stopifnot(length(incidental.x) == length(incidental.weight))
	weight <- new.Weight(focus.weight = focus.weight, incidental.weight = incidental.weight)
	MLE(object = object, x = focus.x, incidental.x = incidental.x, weight = weight, ...)
}
weighted.loglikFun <- function(object, focus.x, incidental.x, focus.weight, incidental.weight, ...)
{# returns the weighted loglikelihood function; user-friendly version of weighted.likFun

	assert.is(is(object, "Family"))
	assert.are(list(focus.x, incidental.x, focus.weight, incidental.weight), "numeric")
	stopifnot(length(focus.x) == length(focus.weight))
	stopifnot(length(incidental.x) == length(incidental.weight))
	focus.fun <- likFun(object = object, x = focus.x, ...)
	weight <- new.Weight(focus.weight = focus.weight, incidental.weight = incidental.weight)
	weighted.likFun(object = object, incidental.x = incidental.x, weight = weight)
}

weighted.likFun <- function(object, incidental.x, weight)
{
	assert.is(object, "likFun")
	assert.is(incidental.x, "numeric")
	assert.is(weight, c("numeric", "Weight"))
	if(is(weight, "numeric") && length(weight) == 1 && length(incidental.x) > 0)
		weight <- rep(weight, length(incidental.x))
	num.incidental.w <- if(is(weight, "numeric"))
		length(weight)
	else if(is(weight, "Weight"))
		length(weight@incidental.weight)
	else
		stop("bad class of weight")
	stopifnot(length(incidental.x) == num.incidental.w)
	stopifnot(length(object@base) == 1 && object@base > 0)
	if(is.nothing(incidental.x) && is.nothing(weight))
	{
		stopifnot(is.nothing(weight))
		object
	}
	else if(!is.nothing(incidental.x) && !is.nothing(weight))
	{
		x.lik.fun <- as(object, "function")
		incidental.x.lik.fun <- function(...)
		{
			nam <- names(object@arglis)
			vec <- sapply(incidental.x, function(wx)
			{
				arglis <- object@arglis[nam != "x"]
				arglis$x <- wx
				fun <- do.call(likFun, arglis)
				sca <- fun(...)
				if(!is.numeric(sca) || length(sca) != 1)
				{ message("sca problem"); browser()}
				sca
			})
			if(!is.numeric(vec))
			{ message("!is.numeric(vec)"); browser()}
			vec
		}
		if(!is(weight, "Weight"))
			weight <- new.Weight(incidental.weight = weight)
		object@weight <- weight
		if(is.nothing(weight)) stop("bad weight")
		object@.Data <- function(...)
		{
			x.lik <- x.lik.fun(...)
			incidental.x.lik <- incidental.x.lik.fun(...)
			assert.is(incidental.x.lik, "numeric")
			assert.is(weight, "Weight")
			stopifnot(length(incidental.x.lik) == length(weight@incidental.weight))
			join <- function(x.part, incidental.x.part)
			{
				stopifnot(length(x.part) == 1)
				vec <- c(x.part, incidental.x.part)
				assert.is(vec, "numeric")
				vec
			}
			normalized.w <- join(x.part = weight@focus.weight, incidental.x.part = weight@incidental.weight)
			combined.lik <- join(x.part = x.lik, incidental.x.part = incidental.x.lik)
			stopifnot(length(normalized.w) == length(combined.lik))
			scalar(sum(normalized.w * combined.lik))
		}
		object
	}
	else
		stop("weighting arguments are bad")
} # end weighted.likFun

hindsight.type <- function(){"hindsight maximum likelihood"}
setMethod("new.log", representation(old.log = "lik.ratio", base = "numeric", old.base = "missing", new.base = "missing"), function(old.log, base, old.base, new.base)
{
	object <- old.log
	object@.Data <- new.log(old.log = as(object, "numeric"), base = base, old.base = object@base)
	object@base <- Scalar(base)
	object
})
codelength <- function(object, hindsight = numeric(0), base = 2, old.base, favoring.null = TRUE)
{
	stopifnot(length(base) == 1 && base > 0)
	if(is(object, "numeric"))
	{
		if(is.nothing(hindsight))
		{
			lr <- if(missing(old.base))
			{
				assert.is(object, "lik.ratio")
				lik.ratio(object, base = base)
			}
			else
				new.log(old.log = object, base = base, old.base = old.base)
			fac <- if(favoring.null) -1 else 1
			vec <- as(fac * lr, "numeric")
			names(vec) <- names(object)
			vec
		}
		else if(is(hindsight, "lik.ratio") && hindsight@type == hindsight.type() && sameNames(object, hindsight) && missing(old.base))
			regret(object = object, hindsight = hindsight, base = base)
		else
			stop("arg err")
	}
	else if(is(object, "lik.ratios")) # added 111022
	{
		get.cl <- function(elem)
		{
			codelength(object = elem, hindsight = hindsight, base = base, favoring.null = favoring.null)
		}
		lapply(object, get.cl)
	}
	else
		stop("bad object")
}
regret <-	function(object, hindsight, base = 2)
{
	assert.is(object, "lik.ratio")
	assert.is(hindsight, "lik.ratio")
	stopifnot(length(base) == 1 && base > 0)
	stopifnot(hindsight@type == hindsight.type())
	stopifnot(sameNames(object, hindsight))
	get.cl <- function(x){codelength(object = x, base = base)}
	get.cl(object) - get.cl(hindsight)
}
single.model.redundancy <-	function(object, max.term = Inf, ...)
{
	if(is(object, "lik.ratio"))
	{
		cl <- codelength(object = object, ...)
		cl <- ifelse(cl <= max.term, cl, max.term)
		scalar(sum(cl))
	}
	else if(is(object, "list"))
	{
		object <- lik.ratios(object)
		re <- sapply(object, single.model.redundancy, max.term = max.term, ...)
		names(re) <- type(object)
		re
	}
	else
		stop("cannot compute single.model.redundancy from the arguments given")
}
redundancy <- function(object1, object2, ...)
{
	assert.is(object1, c("lik.ratio", "list"))
	if(missing(object2))
		single.model.redundancy(object = object1, max.term = 0, ...) # min(a, b) == min(a - b, 0) + b
	else if(are(list(object1, object2), "lik.ratio") && sameNames(object1, object2))
	{
		get.cl <- function(object){codelength(object = object, ...)}
		cl1 <- get.cl(object1)
		cl2 <- get.cl(object2)
		stopifnot(sameNames(cl1, cl2))
		nam <- names(cl1)
		vec <- pmin(cl1, cl2)
		names(vec) <- nam
		stopifnot(sameNames(object1, vec))
		scalar(sum(vec))
	}
	else if(is(object1, "list"))
	{
		object <- lik.ratios(object)
		typ <- type(object1)
		mre <- sapply(1:length(typ), function(i)
		{
			object2.i <- if(is(object2, "lik.ratios") && sameNames(object1, object2) && all(type(object1) == type(object2)))
				object2[[i]]
			else if(is(object2, "lik.ratios") && length(object2) == 1)
				object2[[1]]
			else if(is(object2, "lik.ratio"))
				object2
			else
				stop("cannot compute the minimum codelength from the arguments given")
			redundancy(object1 = object1[[i]], object2 = object2.i, ...)
		})
		names(mre) <- typ
		mre
	}
	else
		stop("bad redundancy args")
} # end redundancy

new.min.redundancy <- function(object, lr, family, unknown.param.1, unknown.param.2)
{
	new("min.redundancy", scalar(object), lr = lr, family = family, unknown.param.1 = unknown.param.1, unknown.param.2 = unknown.param.2)
}
min.redundancy <- function(x, object2, fam, optimize.args = list(), verbose = FALSE, need.estimate = FALSE, scalar.density = FALSE, ...) # family.fun.arg.name,
{
	assert.is(fam, "Family")
	stopifnot(length(fam@unknown.param.name) == 1)
	assert.is(object2, c("list", "univariateDistribution"))
	if(is(object2, "list"))
		stopifnot(all(names(object2) == fam@unknown.param.name))
	assert.is(optimize.args, "list")
	if(is(x, "xprnSetObject") || (!scalar.density && is(x, "numeric")))
	{
		get.MLE <- function(x){MLE(object = fam, x = x)}
		mle <- if(is(x, "xprnSetObject"))
			get.MLE(x = x)
		else if(!scalar.density && is(x, "numeric"))
			vectorize(get.MLE)(x)
		else
			stop("be more careful")
		optimize.args <- c(optimize.args, list(lower = min(mle), upper = max(mle)))
	}
	get.unknown.param <- function(unknown.param.value)
	{
		unknown.param <- list(unknown.param.value)
		names(unknown.param) <- fam@unknown.param.name
		unknown.param
	}
	get.lr <- function(unknown.param.value, verbose)
	{
		unknown.param <- get.unknown.param(unknown.param.value)
		lik.ratio(x = x, object1 = unknown.param, object2 = object2, fam = fam, verbose = verbose, need.estimate = need.estimate, scalar.density = scalar.density, ...) # "scalar.density" added 100911
	}
	get.red <- function(unknown.param.value, verbose = FALSE)
	{
		lr <- get.lr(unknown.param.value, verbose = verbose)
		redundancy(lr)
	}
	opt <- do.call(optimize, c(list(f = get.red, maximum = FALSE), optimize.args))
	unknown.param.value <- opt$minimum
	unknown.param <- get.unknown.param(unknown.param.value)
	lr <- get.lr(unknown.param.value, verbose = verbose)
	if(verbose)
	{
		message("min.redundancy unknown.param: ", unknown.param)
		browser()
	}
	new.min.redundancy(object = opt$objective, lr = lr, family = fam, unknown.param.1 = unknown.param, unknown.param.2 = object2)
}

new.lik.ratio <- function(object, type, base, estimate = numeric(0), unknown.param.set = list())
{
	new("lik.ratio", object, type = type, base = Scalar(base), unknown.param.set = unknown.param.set, estimate = estimate)
}
#XXX|: to MDL.R: nml.type <- "normalized maximum likelihood"
lik.ratio <- function(x, object1, object2, base = exp(1), old.base, type, FUN = NULL, x.converges, fam = NULL, leave1out = FALSE, verbose = FALSE, need.estimate = FALSE, ...)
{
	get.lik.ratio <- function(vec, type, estimate = numeric(0), condition = 0, ...)
	{
		if(need.estimate && is.nothing(estimate))
		{	message(condition, " will not return any estimate"); browser()}
		lr <- new.lik.ratio(object = vec, type = type, base = base, estimate = estimate, ...)
		if(need.estimate)
		{
			if(is.nothing(lr@estimate))
			{	message( "condition ", condition, " failed to produce any estimate"); browser()}
			else if(any(is.na(lr@estimate)))
			{	message( "condition ", condition, " failed produced NA estimate"); browser()}
		}
		if(verbose) {message("lr:"); print(stats(lr))}
		lr
	}
	lr <- if(missing(object1) && !missing(object2) && is(object2, "list") && missing(old.base) && missing(type) && is.null(FUN) && missing(x.converges) && is(fam, "Family"))
	{
		condition <- 1
		if(verbose)
			message("condition ", condition)
		get.mr <- function(y){min.redundancy(x = y, object2 = object2, base = base, fam = fam, verbose = verbose, need.estimate = need.estimate, ...)}
		if(leave1out)
		{
			if(verbose)
				message("leaving one out")
			assert.is(x, "xprnSetObject")
			nam <- names(x)
			estimate <- lr <- as.numeric(rep(NA, length(nam)))
			for(i in 1:length(nam))
			{
				mr <- get.mr(x[-i, ])
				stat.i <- if(is(x, "xprnSet"))
					nstatistic(object = fam, x = exprs(x)[i, ])
				else if(is(x, "xprnSetPair"))
					nstatistic(object = fam, x = exprs(x@x)[i, ], y = exprs(x@y)[i, ])
				else
					stop("no way")
				lik.r <- lik.ratio(x = stat.i, object1 = mr@unknown.param.1, object2 = mr@unknown.param.2, base = base, fam = fam, verbose = verbose, need.estimate = need.estimate)
				lr[i] <- scalar(lik.r)
				if(length(lik.r@estimate) != 1)
				{ message("wrong lik.r@estimate length:"); print(lik.r@estimate); browser()}
				estimate[i] <- scalar(lik.r@estimate)
			}
			if(any(is.na(lr)))
			{ message("missing likelihood ratio"); browser()}
			names(lr) <- nam
			get.lik.ratio(vec = lr, type = "likelihood", estimate = estimate, condition = condition)
		} # end leave1out
		else
			get.mr(y = x)@lr
	}
	else if(missing(object1) && missing(object2) && is.null(FUN) && missing(x.converges) && length(list(...)) == 0 && is.null(fam) && !leave1out)
	{
		condition <- 2
		if(verbose)
			message("condition ", condition)
		if(is(x, "numeric") && !is(x, "lik.ratio") && !missing(old.base) && !missing(type))
			get.lik.ratio(vec = new.log(old.log = x, base = base, old.base = old.base), type = type, condition = condition)
		else if(is(x, "lik.ratio") && missing(old.base) && missing(type))
			new.log(old.log = x, base = base)
		else
		{ message("missing object1 and object2 but bad lik.ratio args"); browser()}
	}
	else if(!missing(x) && is(x, "xprnSetObject") && missing(x.converges) && missing(old.base) && missing(type) && !leave1out)
	{
		condition <- 3
		if(verbose)
			message("condition ", condition)
		nam <- featureNames(x)
		estimate <- vec <- as.numeric(rep(NA, length(nam)))
		type <- unknown.param.set <- NULL
		get.fam <- function()
		{
			fa <- if(is(object1, "familyContaining") && is.null(fam))
				family(object1)
			else if(!is(object1, "familyContaining") && is(fam, "Family"))
				fam
			else
				stop("cannot verify degrees of freedom for expression data using the arguments given")
			fa
		}
		get.df <- function(){get.fam()@known.param$df}
		for(i in 1:length(vec))
		{
			get.row <- function(es)
			{
				assert.is(es, "xprnSet")
				ro <- as(exprs(es)[i, ], "numeric") # corrected 21 October 2011
				if(is(es, "XprnSet"))
					ro <- logb(ro)
				ro
			}
			if(is(x, "xprnSet"))
			{
				x.i <- get.row(es = x)
				assert.is(x.i, "numeric")
				stat <- if(sum(is.finite(x.i)) - 1 == get.df())
				{
					nstatistic(object = get.fam(), x = x.i)
				}
				else
				{	message("not implemented for get.fam()@Distr.fun and xprnSet ", class(x)); print(get.fam()@Distr.fun); browser()}
				assert.is(stat, "statistic")
			}
			else if(is(x, "xprnSetPair"))
			{
#				x@x <- get.es(es = x@x)
#				x@y <- get.es(es = x@y)
				x.i <- get.row(es = x@x)
				y.i <- get.row(es = x@y)
				assert.are(list(x.i, y.i), "numeric")
				stat <- if(sum(is.finite(c(x.i, y.i))) - 2 == get.df())
				{
					nstatistic(object = get.fam(), x = x.i, y = y.i)
				}
				else
				{	message("not implemented for get.fam()@Distr.fun and xprnSetPair ", class(x)); print(get.fam()@Distr.fun); browser()}
				assert.is(stat, "statistic")
			}
			else
				stop("argument failure")
			if(verbose)
			{
				message("computed stat")
				print(class(object1))
				print(class(object2))
			}
			stopifnot(length(stat) == 1)
			sca <- lik.ratio(x = stat, object1 = object1, object2 = object2, base = base, FUN = FUN, fam = fam, verbose = verbose, need.estimate = need.estimate, ...) # exp(1), x.converges = FALSE,  before 100802
#if(is(object1, "familyContaining") && is.null(fam))
#			else
#				stop("")
			stopifnot(length(sca) == 1)
			vec[i] <- sca
			if(length(sca@estimate) != 1)
			{ message("wrong estimate length:"); print(sca@estimate); browser()}
			estimate[i] <- sca@estimate
			type <- sca@type
			unknown.param.set <- sca@unknown.param.set
		} # end for(i in 1:length(vec))
		stopifnot(!is.null(type) && !is.null(unknown.param.set))
		ok <- length(vec) == length(nam)
		if(!ok)
		{ message("bad vec?"); browser()}
		if(!all(is.finite(vec)))
			warning("non-finite likelihood ratios")
		names(vec) <- nam
		get.lik.ratio(vec = vec, type = type, unknown.param.set = unknown.param.set, estimate = estimate, condition = condition)
	}
	else if(!missing(x) && is(x, "numeric") && are(list(object1, object2), c("univariateDistribution", "function", "list")) && is.null(FUN) && missing(x.converges) && missing(old.base) && missing(type) && !leave1out)
	{
		condition <- 4
		if(verbose)
			message("condition ", condition)
		if(is(object1, "list") || is(object2, "list"))
		{
			if(is.null(fam))
				stop("fam is missing from the argument list")
			else
			{
				assert.is(fam, "Family")
				prepare <- function(object)
				{
					if(is(object, "list"))
						Distr(object = fam, unknown.param = object)
					else
						object
				}
				object1 <- prepare(object1)
				object2 <- prepare(object2)
			}
		}
		else
			stopifnot(is.null(fam))
		assert.are(list(object1, object2), c("univariateDistribution", "function"))
		pdens <- function(object){prob.density(object = object, x = x, base = base, ...)}
		vec <- pdens(object = object1) - pdens(object = object2)
		estimate <- if(is(object1, "list") && length(object1) == 1 && is.numeric(object1[[1]]))
			object1[[1]]
		else if(is(object1, "univariateDistribution") && "ncp" %in% slotNames(object1))
			object1@ncp
		else if(need.estimate)
		{ message("the needed estimate is not available"); browser()}
		else
			numeric(0)
		if(verbose)
		{
			message("x:"); print(stats(x))
			message("vec:"); print(stats(vec))
		}
		get.lik.ratio(vec = vec, type = "likelihood", estimate = estimate, condition = condition)
	}
	else if((is(object2, "list") || is(object2, "numeric")) && missing(old.base) && missing(type) && is.null(fam) && !leave1out)
	{
		if(missing(x) && is(object1, "DensityFunNMLs"))
		{
			condition <- 4.7
			if(verbose)
				message("condition ", condition)
			if(missing(x.converges))
				x.converges <- FALSE
			lis <- lapply(1:length(object1), function(i)
			{
				lik.ratio(x = new.statistic(object1@reduced.x)[i], object1 = object1[[i]], object2 = object2, leave1out = leave1out, fam = fam, ...)
			})
			names(lis) <- names(object1)
			vec <- slots(lis)
			estimate <- slots(lis, name = "estimate")
			estimate.ok <- is(estimate, "numeric") && length(estimate) %in% c(0, 1, length(vec))
			if(!estimate.ok)
			{ message("estimate is not suitable for the target class"); browser()}
			get.lik.ratio(vec = vec, type = nml.type, unknown.param.set = if(is(fam, "finiteFamily")) stop("finiteFamily pending implementation") else list(), estimate = estimate, condition = condition)
		} # end condition 4.7
		else if(!missing(x) && is(x, "numeric") && is(object1, "familyContaining"))
		{
			condition <- 5
			if(verbose)
				message("condition ", condition)
			if(missing(x.converges))
				x.converges <- FALSE
			fam <- family(object1) #if(is(object1, "Family"))
			assert.is(fam, "Family")
			if(!is(x, "statistic"))
				x <- nstatistic(fam, x = x, x.converges = x.converges)
			if(length(x) > 1)
				warning("Does x require reduction?")
			lf <- likFun(object = fam, x = x, base = base)
			val0 <- if(is(object2, "list"))
			{
				if(length(object2) == 1)
					object2[[1]]
				else
					stop("bad object2")
			}
			else if(is(object2, "numeric"))
				object2
			else
				stop("cannot get val0")
			assert.is(val0, "numeric")
			stopifnot(length(val0) == 1)
			basic.args.ok <- ((is(object2, "list") && sameSet(names(object2), fam@unknown.param.name)) || is(object2, "numeric"))
			if(!basic.args.ok)
			{ message("bad lik.ratio args 0"); browser()}
			estimate <- numeric(0)
			numer <- if(is(object1, "DensityFunNML") && is.null(FUN))
			{
				nml <- object1(x = x)
				logb(nml, base = base)
			}
			else if(is(object1, "finiteFamily") && is(FUN, "function") && length(object1@unknown.param.set) == 1)
			{
				param.space <- object1@unknown.param.set[[1]]
				assert.is(param.space, "numeric")
				FUN(sapply(param.space, lf))
			}
			else if(is(object1, "extendedFamily") && is.null(FUN))
			{
				mle.and.ml <- MLE(object = object1, base = base, return.MLE = logical(0), log = TRUE, x = x)
				okay <- !is.null(names(mle.and.ml)) && (is.numeric(mle.and.ml) && all(c("mle", "ml") %in% names(mle.and.ml)))
				estimate <- try(mle.and.ml["mle"])
				if(!okay || is.err(estimate))
				{ message("mle.and.ml not okay:"); print(mle.and.ml); browser()}
				mle.and.ml["ml"]
			}
			else
			{ message("bad lik.ratio args 1"); browser()}
			vec <- scalar(numer - lf(val0))
			type <- if(is(object1, "finiteFamily") && all(sapply(object1@unknown.param.set, length) == 1))
				"likelihood"
			else if(is(object1, "DensityFunNML"))
				nml.type
			else if(identical(FUN, mean))
				"Bayes factor"
			else if(identical(FUN, max) || is.null(FUN))
				hindsight.type()
			else
				"unrecognized type"
			get.lik.ratio(vec = vec, type = type, unknown.param.set = if(is(fam, "finiteFamily")) fam@unknown.param.set else list(), estimate = estimate, condition = condition)
		} # end condition 5
		else
			stop("subconditions do not hold")
	} # end conditions 4.7 and 5
	else
	{ message("bad lik.ratio args 2"); browser()}
	assert.is(lr, "lik.ratio")
	lr
} # end lik.ratio
lik.ratio.ibe <- function(x, ...)
{
	assert.is(x, "xprnSetObject")
	lik.ratio(x = x, object1 = extendedFamily.absTd(df = x), object2 = list(ncp = 0), ...)
}
Lik.ratio <- function(...)
{
	lr <- lik.ratio(..., base = exp(1))
#	new.lik.ratio(object = exp(lr), type = type, base = base)
	lr@.Data <- exp(lr)
	lr@base <- Scalar(numeric(0))
	lr
}
lik.ratio.Td <- function(x, df, ncp, ncp.set, family.fun = Family.Td, finite.family.fun = finiteFamily.Td, ncp0 = 0, mixture = FALSE, verbose = FALSE, save.time = FALSE, individual = save.time, ...) # , P0 = numeric(0)
{
	if(verbose)
		message("lik.ratio.Td starting on ", date())
	if(missing(df))
	{
		assert.is(x, "xprnSetObject")
		df <- Df(x)
	}
	assert.is(df, "numeric")
	get.arglis <- function(...)
	{
		arglis <- list(...)
		if(mixture)
			arglis <- c(arglis, list(ncp0 = ncp0))
		arglis
	}
	if(missing(ncp) && missing(ncp.set))
	{
		arglis.proper <- list(...)
		recall <- function(save.time, ...)
		{
			do.call("lik.ratio.Td", c(arglis.proper, list(x = x, df = df, family.fun = family.fun, finite.family.fun = finite.family.fun, ncp0 = ncp0, mixture = mixture, verbose = verbose, save.time = save.time, ...)))
		}
		get_statistic <- function()
		{
			sta <- if(is(x, "numeric"))
				x
			else if(is(x, "xprnSetObject"))
				nstatistic(x)
			else
				stop("cannot get fast ncp")
			stopifnot(length(sta) == Length(x))
			sta
		}
		get.fam <- function()
		{
			arglis <- get.arglis(df = df)
			fam <- try(do.call(family.fun, arglis))
			if(is.err(fam))
			{ message("first family problem"); browser()}
			fam
		}
		if(individual && !save.time) # added 111022
		{
			sta <- get_statistic()
			fam <- get.fam()
			ncp <- sapply(sta, function(sca){as(MLE(fam, x = sca), "numeric")})
			lr <- try(recall(save.time = FALSE, ncp = as(ncp, "numeric")))
			lr@estimate <- ncp
			ok <- is(lr, "lik.ratio") && length(lr) == length(lr@estimate)
			if(!ok)
			{ message("bad lr"); browser()}
			names(lr@estimate) <- names(lr)
		}
		else if(individual && save.time) # added 111022
		{
			if(verbose) message("lik.ratio.Td condition A1")
			ncp <- get_statistic()
			lr <- recall(save.time = FALSE, ncp = ncp)
		}
		else if(!individual && !save.time)
		{
			if(verbose) message("lik.ratio.Td condition A2")
			fam <- get.fam()
			lr <- try(lik.ratio(x = x, object2 = list(ncp = ncp0), fam = fam, verbose = verbose, ...), silent = FALSE)
			if(is.err(lr)) # added 111022
			{
				if(is(fam, "Family") && !mixture)
				{
					warning("switched to individual = TRUE")
					lr <- try(recall(save.time = FALSE, individual = TRUE), silent = FALSE)
				}
				else
				{
					warning("switched to save.time = TRUE")
					lr <- try(recall(save.time = TRUE), silent = FALSE)
				}
			}
		}
		else
			stop("bad individual & save.time")
		if(is.err(lr))
		{ message("lr error"); browser()}
		lr
	}
	else if(!individual && !save.time && missing(ncp))
	{
		if(verbose) message("lik.ratio.Td condition B")
		fam <- try(finite.family.fun(df = df, ncp.set = ncp.set))
		if(is.err(fam))
		{ message("family problem"); browser()}
		lik.ratio(x = x, object1 = fam, object2 = list(ncp = ncp0), FUN = max, verbose = verbose, ...)
	}
	else if(!individual && !save.time && missing(ncp.set))
	{
		if(verbose)
		{
			message("lik.ratio.Td condition C")
		  message("ncp: ")
		  print(ncp)
		  message("")
		}
		fam <- family.fun(df = df)
		get.lr <- function(y, ncp1)
		{
			stopifnot(length(ncp1) == 1)
			if(verbose)
				print(list(x = y, object1 = list(ncp = ncp1), object2 = list(ncp = ncp0), fam = fam, verbose = verbose, ...))
			lik.ratio(x = y, object1 = list(ncp = ncp1), object2 = list(ncp = ncp0), fam = fam, verbose = verbose, ...)
		}
		if(length(ncp) == 1)
			get.lr(y = x, ncp1 = ncp)
		else if(length(ncp) == Length(x)) # added 21 October 2011
		{
			get.y <- function(i)
			{
				if(is(x, "xprnSet"))
					x[i, ]
				else if(is(x, "numeric"))
					x[i]
				else
					stop("cannot get y")
			}
			vec <- estimate <- as.numeric(rep(NA, length(ncp)))
			for(i in c(1:length(ncp)))
			{
				lr <- get.lr(y = get.y(i), ncp1 = ncp[i])
				if(i == 1)
				{
					base <- lr@base
					type <- lr@type
				}
				estimate[i] <- lr@estimate
				vec[i] <- as(lr, "numeric")
			}
			stopifnot(length(vec) == length(estimate))
			names(vec) <- names(estimate) <- names(x)
			new.lik.ratio(object = vec, type = type, base = base, estimate = estimate)
		}
		else
		{	message("bad 21 October 2011"); browser()}
	}
	else
		stop("bad lik.ratio.Td")
} # end lik.ratio.Td
lik.ratio.absTd <- function(...){lik.ratio.Td(..., family.fun = extendedFamily.absTd, finite.family.fun = finiteFamily.absTd)} # family.fun = Family.absTd before 100806
lik.ratio.absTdMixture <- function(x, estimate, family.fun = extendedFamily.absTd, df, ncp0 = 0, base = exp(1), ...)
{# {lik.ratio.Td(..., family.fun = extendedFamily.absTdMixture, finite.family.fun = finiteFamily.absTdMixture, mixture = TRUE)}

	if(missing(df))
	{
		if(is(x, "xprnSetObject"))
			df <- Df(x)
		else
			stop("cannot get df from x")
	}
	stopifnot(is(df, "numeric") && length(df) == 1)
	fam <- family.fun(df = df)
	if(missing(estimate))
		stop("estimate is missing; in principle, it could be extracted from x, but that has not been implemented") # object <- posteriorP0()
	if(is(estimate, "posteriorP0"))
		estimate <- estimate@estimate
	assert.is(estimate, "nmle")#XXX|:changing from mle to nmle
	if(is(x, "xprnSetObject"))
		x <- nstatistic(object = fam, x = x)
	assert.is(x, "numeric")
	get.object <- function(ncp)
	{
		assert.is(ncp, "numeric")
		stopifnot(length(ncp) == 1)
#		dis <- Distr(object = fam, unknown.param = list(ncp = ncp))
#		stopifnot(is(dis, "univariateDistribution"))
#		dis
		list(ncp = ncp)
	}
	ncp1 <- estimate@unknown.param$ncp
	lr <- sapply(x, function(sca){lik.ratio(x = scalar(sca), object1 = get.object(ncp = ncp1), object2 = get.object(ncp = ncp0), fam = fam, ...)})
	stopifnot(length(lr) == length(x))
	lr <- new.lik.ratio(lr, type = "Bayes factor estimate", base = base)
	names(lr) <- names(x)
	lr
}

lik.ratios <- new.lik.ratios <- function(object, base)
{
	if(missing(base) && is(object, "lik.ratios"))
		object
	else
	{
		if(missing(base))
			base <- 2
		if(!missing(object) && length(object) > 0)
		{
			assert.is(object, "list")
			lrs <- lapply(object, lik.ratio, base = base)
			new("lik.ratios", lrs)
		}
		else
			stop("new.lik.ratios arg error")
	}
}
base <- function(object)
{
	if(is(object, "lik.ratios"))
		base(object[[1]])
	else
		object@base
}
type <- function(object)
{
	if(is(object, "lik.ratios"))
		sapply(object, slot, name = "type")
	else
		object@type
}

new.mle <- function(object = numeric(0), unknown.param = list())
{#XXX|:changing from mle to nmle
	if(is(object, "nmle") && (length(unknown.param) == 0 || identical(unknown.param, object@unknown.param)))
		object
	else if(is(object, "numeric"))#XXX|:changing from mle to nmle
		new("nmle", object, unknown.param = unknown.param)
	else
		stop("cannot create nmle")
}
mle <- function(object, call.browser = FALSE, par = numeric(0), optim.args = list(), ...) # save.time = !is(object, "extendedFamily") || length(optim.args) > 0
{
	if(call.browser)
	{ message("list(...): "); print(list(...)); browser()}
	assert.is(object, "Family")
	if(is.nothing(par))
	{
		assert.is(object, "extendedFamily")
		cha <- object@mle.fun(...)
		new.mle(cha)
	}
	else
	{
		assert.is(par, "numeric")
		assert.is(names(par), "character")
		lfun <- likFun(object = object, ...)
		named.list <- function(vec)
		{
			assert.is(vec, "numeric")
			stopifnot(length(vec) == length(par))
			names(vec) <- names(par)
			as.list(vec)
		}
		fn <- function(vec)
		{
			unknown.param <- named.list(vec)
			if("P0" %in% names(unknown.param) && unknown.param$P0 > 1)
			{
				# warning("    unknown.param$P0 of ", round(unknown.param$P0, 2), " reset to 1 on ", date())
				unknown.param$P0 <- 1
			}
			lik <- try(do.call(lfun, unknown.param))
			assert.is(lik, "numeric")
			stopifnot(length(lik) == 1)
			lik
		} #(P0 = object[1], ncp = object[2])
		get.opt <- function(pa) # pa = c
		{
			do.call(optim, c(list(par = pa, fn = fn, control = list(fnscale = -1)), optim.args))
		}
		opt <- get.opt(pa = par)
		if(opt$convergence != 0 && "P0" %in% names(par) && "ncp" %in% names(par))
		{
			pa <- par
			pa["P0"] <- 1 # (par["P0"] + 9) / 10
			pa["ncp"] <- 0 # par["ncp"] / 2 # 2 * par["ncp"]
			opt <- get.opt(pa = pa)
		}
		if(opt$convergence != 0)
		{ message("failure to converge: ", opt$convergence); browser()}
		unknown.param <- named.list(opt$par)
		ml <- opt$value
		new.mle(c(ml = ml, mle = opt$par[1]), unknown.param = unknown.param)
	}
} # end mle

MLE <- function(object, lower, upper, lower.unknown.param, upper.unknown.param, return.MLE = TRUE, log = FALSE, base = exp(1), x, truncation.tol = numeric(0), boundary.tol = numeric(0), incidental.x = numeric(0), weight = numeric(0), excuse.infinities = TRUE, browser.range = numeric(0), ...) # excuse.infinities = logical(0)
{
	# print(browser.range); message("browser.range check A"); browser()
	ml.return.value <- function(ml, Base)
	{
		return.value <- if(log)
			ml
		else
			Base ^ ml
		ok <- (log && return.value == -Inf) || is.finite(return.value)
		if(!ok)
		{ message("return.value error"); browser()}
		return.value
	}
	finite.parameter.set <- try(is(object, "finiteFamily"))
	if(is.err(finite.parameter.set))
	{ message("cannot determine whether parameter set is finite"); browser()}
	get.weighted.likFun <- function(object)
	{
		assert.is(object, "likFun")
		fun <- weighted.likFun(object = object, incidental.x = incidental.x, weight = weight)
		assert.is(fun, "likFun")
		fun
	}
	to.return <- if(!missing(x) && is(x, "xprnSetObject") && !missing(object) && is(object, "Family") && missing(lower) && missing(upper) && missing(lower.unknown.param) && missing(upper.unknown.param))
	{
		get_statistic <- function(...){nstatistic(object = object, ...)}
		FUN <- if(is(x, "xprnSet"))
			function(x, na.rm)
			{
				get_statistic(x = x)
			}
		else if(is(x, "xprnSetPair"))
			function(x, y, na.rm)
			{
				get_statistic(x = x, y = y)
			}
		else
			stop("that cannot be")
		sta <- stat(x = x, FUN = FUN)
		vec <- sapply(sta, function(st){MLE(x = st, object = object, return.MLE = return.MLE, log = log, base = base, truncation.tol = truncation.tol, incidental.x = incidental.x, weight = weight, browser.range = browser.range, ...)})
		names(vec) <- names(x)
		vec
	}
	else if(finite.parameter.set && missing(lower) && missing(upper) && missing(lower.unknown.param) && missing(upper.unknown.param))
	{
		f <- get.weighted.likFun(likFun(object = object, base = base, x = x))
		domain <- if(length(object@unknown.param.set) == 1)
			object@unknown.param.set[[1]]
		else
		{	message("cannot determine the domain of the likelihood function"); browser()}
		opt <- finite.optimize(f = f, domain = domain, maximum = TRUE, ...)
		mle <- opt$maximum
		ml <- ml.return.value(ml = opt$objective, Base = f@base)
		if(is.nothing(return.MLE))
			c(mle = mle, ml = ml)
		else if(return.MLE)
			mle
		else
			ml
	}
	else if(is(object, "likFun") && missing(x) && !missing(lower) && !missing(upper))
	{
		stopifnot(length(object@unknown.param.name) == 1)
		end.point <- function(bound)
		{
			if(is(bound, "function"))
				bound <- bound(object@x)
			else if(!is.finite(bound))
			{ message("non-finite bound"); browser()}
			bound.for.optimize <- if(is(bound, "numeric") && length(bound) == 1)
				bound
			else
				stop("bad bound in end.point")
			stopifnot(is.finite(bound.for.optimize))
			bound.for.optimize
		}
		if("return.mle" %in% names(list(...)) || "return.MLE" %in% names(list(...)))
		{ message("invalid argument:"); print(names(list(...))); browser()}
		lower.boundary <- end.point(lower)
		upper.boundary <- end.point(upper)
		weighted.object <- get.weighted.likFun(object = object)
		assert.is(weighted.object, "likFun")
		opt <- if(lower.boundary <= upper.boundary)
		{
			try(Optimize(weighted.object, lower = lower.boundary, upper = upper.boundary, maximum = TRUE, lower.reliable = FALSE, boundary.tol = boundary.tol, excuse.infinities = excuse.infinities, ...))
		}
		else
			stop("bad boundaries")
		if(is.err(opt) || !is(opt, "list"))
		{
			message("MLE could not maximize the likelihood function")
			plot.lik <- function(length.out = 100, ...)
			{
				X <- seq(lower.boundary, upper.boundary, length.out = length.out)
				Y <- sapply(X, weighted.object)
				plot(x = X, y = Y, xlab = "parameter value", ylab = "weighted.object", ...)
			}
			plot.lik()
			browser()
		}
		mle <- opt$maximum
		if(mle == lower.boundary)
			warning(paste("mle at lower.boundary of", lower.boundary))
		if(mle == upper.boundary)
			warning(paste("mle at upper.boundary of", upper.boundary))
		Lik <- exp(opt$objective)
		if(!is.nothing(browser.range) && Lik >= browser.range[1] && Lik <= browser.range[2])
		{ message("exp(opt$objective) is ", Lik); browser()}
		if(!is.nothing(return.MLE) && return.MLE)
			mle
		else
		{
			ml <- opt$objective
			if(length(truncation.tol) == 1)
			{
				stopifnot(truncation.tol >= 0)
				bound.ok <- function(bound){is(bound, "numeric") && length(bound) == 1}
				if(!bound.ok(lower) || !bound.ok(upper))
				{ message("Is there a problem with bound?"); message()}
				min.mle <- lower.unknown.param + truncation.tol
				max.mle <- upper.unknown.param - truncation.tol
				if(min.mle >= max.mle)
				{ message("impossible bounds"); print(c(min.mle, max.mle)); browser()}
				mle.inside.bounds <- mle >= min.mle && mle <= max.mle
				if(!mle.inside.bounds)
					ml <- -Inf
			}
			else if(length(truncation.tol) > 1)
				stop("bad truncation.tol")
			ml <- ml.return.value(ml = ml, Base = object@base)
			if(is.nothing(return.MLE))
				c(mle = mle, ml = ml)
			else if(!return.MLE)
				ml
			else
				stop("19 not be")
		}
	} # end if(is(object, "likFun") && missing(x) && !missing(lower) && !missing(upper))
	else if(is(object, "extendedFamily") && is.numeric(x) && missing(lower.unknown.param) && missing(upper.unknown.param))
	{
		if(missing(lower))
			lower <- object@lower.unknown.param
		if(missing(upper))
			upper <- object@upper.unknown.param
		object@mle.fun(x = x, return.MLE = return.MLE, log = log, base = base, truncation.tol = truncation.tol, incidental.x = incidental.x, weight = weight, ...)
	}
	else if(is(object, "Family") && is.numeric(x) && missing(lower.unknown.param) && missing(upper.unknown.param))
	{
		lf <- likFun(object = object, base = base, x = x)
		if(missing(lower))
			lower <- object@lower.unknown.param
		if(missing(upper))
			upper <- object@upper.unknown.param
		MLE(object = lf, lower = lower, upper = upper, lower.unknown.param = object@lower.unknown.param, upper.unknown.param = object@upper.unknown.param, return.MLE = return.MLE, log = log, base = base, truncation.tol = truncation.tol, incidental.x = incidental.x, weight = weight, browser.range = browser.range, ...)
	}
	else
		stop("bad MLE args")
	ok <- !is.nothing(return.MLE) || (is.numeric(to.return) && all(c("mle", "ml") %in% names(to.return)))
	if(!ok)
	{ message("data returned at you"); browser()}
	to.return
} # end MLE
extremum.unknown.param <- function(object, name)
{
	extremum <- slot(object, name = name)
	ext <- if(is(object, "finiteFamily") && length(object@unknown.param.set) == 1)
	{
		extreme <- if(name == "upper.unknown.param")
			max
		else if(name == "lower.unknown.param")
			min
		else
			stop("extreme err")
		extreme(object@unknown.param.set[[1]])
	}
	else if(length(extremum) == 1 && is(object, "Family"))
		extremum
	else
		stop("extremum.unknown.param arguments are bad")
	scalar(ext)
}
lower.unknown.param <- function(object){extremum.unknown.param(object, name = "lower.unknown.param")}
upper.unknown.param <- function(object){extremum.unknown.param(object, name = "upper.unknown.param")}
integration.limits <- function(object, p = 1e-3)
{
	if(is(object, "Family"))
	{
		is.ok <- function(x){length(x) == 0 || is.finite(x)}
		stopifnot(is.ok(object@lower.unknown.param) && is.ok(object@upper.unknown.param) && length(p) == 1)
		get.dis <- function(unknown.param.value){Distr(object = object, unknown.param = list(unknown.param.value))}
		get.quantiles <- function(dis)
		{
			qfun <- q(dis)
			c(qfun(p = p), qfun(p = 1 - p))
		}
		quantiles <- sort(c(get.quantiles(get.dis(lower.unknown.param(object))), get.quantiles(get.dis(upper.unknown.param(object)))))
		stopifnot(length(quantiles) == 4)
		lims <- c(lower = quantiles[1], upper = quantiles[4]) # min(quantiles), upper = max(quantiles))
		stopifnot(lims[1] < lims[2])
		lims
	}
	else if(is(object, "DensityFunNML"))
		integration.limits(object = object@family, p = p)
	else
		stop("bad args for integration.limits")
}

Family.Binom <- function(size, lower.prob = 0, upper.prob = 1)
{
	stopifnot(lower.prob < upper.prob)
	Family(object = Binom, unknown.param.name = "prob", known.param = list(size = size), discrete = TRUE, lower.unknown.param = lower.prob, upper.unknown.param = upper.prob)
}
extendedFamily.Binom <- function(size, ...)
{
	object <- Family.Binom(size = size, ...)
	mle.fun <- function(x, ...)
	{
		MLE(object = object, x = x, ...)
	}
	extendedFamily(object = object, mle.fun = mle.fun)
}

Family.Chisq <- function(df, lower.ncp = 0, upper.ncp = Inf, object = Chisq, ncp0, unknown.param.name = "ncp")
{
	stopifnot(lower.ncp >= 0)
	if(upper.ncp <= lower.ncp)
	{ message("Chisq lower.ncp ", lower.ncp, " is not less than upper.ncp ", upper.ncp); browser()}
	known.param <- list(df = df)
	if(!missing(ncp0))
		known.param <- c(known.param, list(ncp0 = ncp0))
	Family(object = object, unknown.param.name = unknown.param.name, known.param = known.param, discrete = FALSE, lower.unknown.param = lower.ncp, upper.unknown.param = upper.ncp)
}
Family.Td <- function(df, lower.ncp = -Inf, upper.ncp = Inf, object = Td, ncp0, unknown.param.name = "ncp")
{
	if(upper.ncp <= lower.ncp)
	{ message("Td lower.ncp ", lower.ncp, " is not less than upper.ncp ", upper.ncp); browser()}
	known.param <- list(df = df)
	if(!missing(ncp0))
		known.param <- c(known.param, list(ncp0 = ncp0))
	Family(object = object, unknown.param.name = unknown.param.name, known.param = known.param, discrete = FALSE, lower.unknown.param = lower.ncp, upper.unknown.param = upper.ncp)
}
Family.absTd <- function(df, lower.ncp = -Inf, upper.ncp = Inf)
{
	Family.Td(df, lower.ncp = -Inf, upper.ncp = Inf, object = absTd)
}
Family.TdMixture <- function(df, lower.ncp = -Inf, upper.ncp = Inf, object = stop("TdMixture not implemented"), ncp0)
{
	if(upper.ncp <= lower.ncp)
	{ message("Family.TdMixture lower.ncp ", lower.ncp, " is not less than upper.ncp ", upper.ncp); browser()}
	known.param <- list(df = df)
	if(!missing(ncp0))
	{
		assert.is(ncp0, "numeric")
		known.param <- c(known.param, list(ncp0 = ncp0))
	}
	Family(object = object, unknown.param.name = c("ncp", "P0"), known.param = known.param, discrete = FALSE, lower.unknown.param = lower.ncp, upper.unknown.param = upper.ncp)
}
Family.TdMixtureGivenP0 <- function(df, lower.ncp = -Inf, upper.ncp = Inf, object = stop("not implemented"), ncp0, P0)
{
	if(upper.ncp <= lower.ncp)
	{ message("Family.TdMixtureGivenP0 lower.ncp ", lower.ncp, " is not less than upper.ncp ", upper.ncp); browser()}
	known.param <- list(df = df, P0 = P0)
	if(!missing(ncp0))
	{
		assert.is(ncp0, "numeric")
		known.param <- c(known.param, list(ncp0 = ncp0))
	}
	Family(object = object, unknown.param.name = c("ncp"), known.param = known.param, discrete = FALSE, lower.unknown.param = lower.ncp, upper.unknown.param = upper.ncp)
}
Family.absTdMixture <- function(df, ncp0 = 0, lower.ncp = -Inf, upper.ncp = Inf) # list(ncp0 = 0)
{
	Family.TdMixture(df, lower.ncp = lower.ncp, upper.ncp = upper.ncp, object = absTdMixture, ncp0 = ncp0)
}
Family.absTdMixtureGivenP0 <- function(df, P0, ncp0 = 0, lower.ncp = -Inf, upper.ncp = Inf) # list(ncp0 = 0)
{
	Family.TdMixtureGivenP0(df, lower.ncp = lower.ncp, upper.ncp = upper.ncp, object = absTdMixture, ncp0 = ncp0, P0 = P0)
}
Family.ChisqMixture <- function(df, ncp0 = 0, lower.ncp = -Inf, upper.ncp = Inf) # list(ncp0 = 0)
{
	Family.TdMixture(df, lower.ncp = lower.ncp, upper.ncp = upper.ncp, object = ChisqMixture, ncp0 = ncp0)
}
finiteFamily.Td <- function(..., ncp.set, family.fun = Family.Td)
{
	assert.is(ncp.set, "numeric")
	new.finiteFamily(family.fun(...), unknown.param.set = list(ncp = ncp.set))
}
finiteFamily.absTd <- function(...){finiteFamily.Td(..., family.fun = Family.absTd)}
finiteFamily.absTdMixture <- function(...){finiteFamily.Td(..., family.fun = Family.absTdMixture)}
Family.Norm <- function(sd, lower.mean = -Inf, upper.mean = Inf, object = Norm, mean0, unknown.param.name = "mean")
{
	if(upper.mean <= lower.mean)
	{ message("Family.Norm lower.mean ", lower.mean, " is not less than upper.mean ", upper.mean); browser()}
	known.param <- list(sd = sd)
	if(!missing(mean0))
		known.param <- c(known.param, list(mean0 = mean0))
	Family(object = object, unknown.param.name = unknown.param.name, known.param = known.param, discrete = FALSE, lower.unknown.param = lower.mean, upper.unknown.param = upper.mean)
}
extendedFamily.Td <- function(df, family.fun = Family.Td, div.mult.factor = numeric(0), lower = 0, browser.range = numeric(0), ...)
{
	if(is.nothing(div.mult.factor))
		div.mult.factor <- 2
	if(is(df, "xprnSetObject"))
		df <- Df(df)
	assert.is(df, "numeric")
	object <- family.fun(df = df, ...)
	lower.ncp <- object@lower.unknown.param
	upper.ncp <- object@upper.unknown.param
	stopifnot(div.mult.factor >= 1)
	mle.fun <- function(x, ...)
	{
		stopifnot(length(x) == 1)
		div <- x / div.mult.factor
		mult <- div.mult.factor * x
		if(is.nothing(lower))
			lower <- max(lower.ncp, min(div, mult))
		upper <- min(upper.ncp, max(div, mult))
		if(upper <= lower)
		{
			if(lower.ncp < lower && is.finite(lower.ncp))
				lower <- lower.ncp
			if(upper.ncp > upper && is.finite(upper.ncp))
				upper <- upper.ncp
		}
		if(lower == 0 && upper == 0)
		{
			lower = -1
			upper = 1
		}
		if(upper < lower) # changed 101129
		{ message("lower ", lower, " is greater than upper ", upper); browser()}
		lower.mle <- lower
		upper.mle <- upper
		MLE(object = object, x = x, lower = lower.mle, upper = upper.mle, browser.range = browser.range, ...)
	}
	extendedFamily(object = object, mle.fun = mle.fun) #, lower.mle = lower.mle, upper.mle = upper.mle)
} # end extendedFamily.Td
extendedFamily.TdMixture <- function(df, family.fun = stop("Family.TdMixture not yet implemented"), lower.ncp = NULL, upper.ncp = NULL, lower.P0 = NULL, upper.P0 = NULL, assumed.P0 = NULL, browser.range = numeric(0), save.time = logical(0), ...) # save.time = logical(0) on 110417; was FALSE
{
	stopifnot(is.null(assumed.P0) || (is.null(lower.P0) || is.null(upper.P0)))
	if(is.null(lower.P0))
	{
		lower.P0 <- if(is.null(assumed.P0)) 0 else assumed.P0
	}
	if(is.null(upper.P0))
	{
		upper.P0 <- if(is.null(assumed.P0)) 1 else assumed.P0
	}
	bounds.set <- !is.null(lower.P0) && !is.null(upper.P0)
	if(!bounds.set)
	{ message("!bounds.set"); browser()}
	if(!is.null(assumed.P0) && !all(assumed.P0 == c(lower.P0, upper.P0)))
	{ message("P0 confusion"); browser()}
	if(is(df, "xprnSetObject"))
		df <- Df(df)
	assert.is(df, "numeric")
	fam <- family.fun(df = df, ...)
	if(is.null(lower.ncp))
		lower.ncp <- 1 # prevents identifiability and numeric problems
	if(is.null(upper.ncp))
		upper.ncp <- fam@upper.unknown.param
	stopifnot(lower.ncp <= upper.ncp)
	stopifnot(lower.ncp >= fam@lower.unknown.param && upper.ncp <= fam@upper.unknown.param)
	fast.mle.fun <- function(x, ...)
	{
		assert.is(x, "numeric")
		stopifnot(length(x) >= 1)
		par <- c(P0 = (upper.P0 + lower.P0) / 2, ncp = mean(x))
		mle(object = fam, x = x, par = par, ...)
	}
	slow.mle.fun <- function(x, return.MLE, base = exp(1), lower = lower.ncp, upper = upper.ncp, verbose, boundary.tol = -Inf, ...)
	{
		if(missing(return.MLE))
			return.MLE = logical(0)
		if(missing(verbose))
			verbose = FALSE
		stopifnot(length(return.MLE) == 0)
		assert.is(x, "numeric")
		stopifnot(length(x) >= 1)
		set.lower <- lower <= fam@lower.unknown.param || !is.finite(lower)
		set.upper <- upper >= fam@upper.unknown.param || !is.finite(upper)
		if(set.lower || set.upper)
		{
			fam1 <- if(identical(family.fun, Family.absTdMixture))
				extendedFamily.absTd(df = df)
			else if(identical(family.fun, Family.ChisqMixture))
				extendedFamily.Chisq(df = df)
			else
			{
				message("cannot compute default lower.ncp or upper.ncp")
				browser()
			}
			if(verbose) message("\n  starting to compute individual ncp estimates for lower or upper on ", date())
			extreme.x <- if(identical(family.fun, Family.absTdMixture) || identical(family.fun, Family.ChisqMixture))
			{
				stopifnot(all(x >= 0))
				range(x)
			}
			else
			{
				message("cannot compute extreme.x")
				browser()
			}
			ncp.vec <- if(length(x) >= 0)
				sapply(if(length(x) > 2) extreme.x else x, fam1@mle.fun)
			else
				stop("bad length(x)")
			stopifnot(length(ncp.vec) %in% c(1, 2))
			if(verbose) message("  finished computing individual ncp estimates for lower or upper on ", date())
			if(set.lower)
				lower <- max((4/5) * min(ncp.vec), lower.ncp)
			if(set.upper)
				upper <- max(lower.ncp, min((5/4) * max(ncp.vec), upper.ncp))
			if(verbose) message("  lower: ", lower, "; upper: ", upper)
		}
		stopifnot(all(is.finite(c(lower = lower, upper = upper))))
		stopifnot(lower <= upper && lower >= fam@lower.unknown.param && upper <= fam@upper.unknown.param)
		if(length(x) == 1 && lower > 0 && upper > 0 && length(ncp.vec) == 1 && ncp.vec >= 0)
		{
			ncp.sca <- Scalar(ncp.vec)
			if(ncp.sca < lower)
			{
				P0.mle <- 1
				ncp.mle <- lower
			}
			else if(ncp.sca >= lower && ncp.sca <= upper)
			{
				P0.mle <- lower.P0
				ncp.mle <- ncp.sca
			}
			else if(ncp.sca > upper)
			{
				P0.mle <- lower.P0
				ncp.mle <- upper
			}
			else
				stop("that did not just happen")
			ml <- as.numeric(NA) # numeric(0) # value can be obtained but is not needed
		}
		else if(length(x) >= 2)
		{
			MLE.given.P0 <- function(P0)
			{
				P0 <- as(P0, "numeric")
				P0.ok <- is.prob(P0)
				if(!P0.ok )
				{ message("!P0.ok; P0: ", P0); browser()}
				fun.of.ncp <- function(ncp)
				{
					dis <- fam@Distr.fun(df = df, P0 = P0, ncp0 = fam@known.param$ncp0, ncp = ncp)
					prob.density(object = dis, x = x, base = base, scalar.density = TRUE)
				}
				lik.fun.of.ncp <- new.likFun(object = fun.of.ncp, Distr.fun = fam@Distr.fun, unknown.param.name = "ncp", known.param = fam@known.param, x = x, base = base, arglis = list())
				MLE(object = lik.fun.of.ncp, lower = lower, upper = upper, browser.range = browser.range, return.MLE = logical(0), boundary.tol = boundary.tol, ...) # returns c(mle = mle, ml = ml)
			}
			get.ml <- function(P0){MLE.given.P0(P0)["ml"]}
			fun.of.P0 <- function(P0)
			{
				fun <- function(p0)
				{
					ml <- get.ml(p0)
					assert.is(ml, "numeric")
					stopifnot(length(ml) == 1)
					ml
				}
				sapply(P0, fun)
			}
			stopifnot(length(lower.P0) == 1 && length(upper.P0) == 1)
			if(lower.P0 == upper.P0)
			{
				P0.mle <- lower.P0
				mle.object <- MLE.given.P0(P0.mle)
				ncp.mle <- mle.object["mle"]
				ml <- mle.object["ml"]
			}
			else
			{
				opt.over.P0 <- try(optimize(fun.of.P0, lower = lower.P0, upper = upper.P0, maximum = TRUE))
				if(is.err(opt.over.P0))
				{ message("opt.over.P0 error"); browser()}
				P0.mle <- opt.over.P0$maximum
				ncp.mle <- MLE.given.P0(P0 = P0.mle)["mle"]
				ml <- opt.over.P0$objective
			}
		} # end if(length(x) >= 2)
		else
			stop("bad length(x) or bad sign")
		names(ml) <- names(ncp.mle) <- NULL
		mle.vec <- c(mle = ncp.mle, ml = ml)
		if(FALSE) # lower.P0 == 0 && lower <= 1e-3 + 1e-4 && P0.mle >= 0.05)
		{ message("I seem to be having trouble finding the maximum likelihood estimate"); browser()}
		new.mle(object = mle.vec, unknown.param = list(ncp = ncp.mle, P0 = P0.mle))
	} # end slow.mle.fun
	mle.fun <- if(is.nothing(save.time))
		function(x, ...)
		{
			mle.obj <- try(slow.mle.fun(x = x, ...))
			if(is.err(mle.obj))
			{
				message("    error from slow.mle.fun; calling fast.mle.fun")
				mle.obj <- fast.mle.fun(x = x, ...)
			}
			mle.obj
		}
	else if(save.time)
	{
		fast.mle.fun
	}
	else
	{
		slow.mle.fun
	}
	extendedFamily(object = fam, mle.fun = mle.fun) #, lower.mle = lower.mle, upper.mle = upper.mle)
} # end extendedFamily.TdMixture

extendedFamily.Norm <- function(sd, family.fun = Family.Norm, null.value = 0, ...)
{
	object <- family.fun(sd = sd, ...)
	mle.fun <- function(x, ...)
	{
		stopifnot(length(x) == 1)
		lower <- x - sd
		upper <- x + sd
		if(upper <= lower)
		{ message("extendedFamily.Norm lower ", lower, " is not less than upper ", upper); browser()}
		lower.mle <- lower
		upper.mle <- upper
		if(!is.nothing(null.value))
		{
			dis <- function(bound){abs(bound - null.value)}
			if(dis(lower.mle) < dis(upper.mle))
				lower.mle <- null.value
			else if(dis(lower.mle) > dis(upper.mle))
				upper.mle <- null.value
			else
				warning("null.value not used by extendedFamily.Norm")
		}
		MLE(object = object, x = x, lower = lower.mle, upper = upper.mle, ...)
	}
	extendedFamily(object = object, mle.fun = mle.fun) #, lower.mle = lower.mle, upper.mle = upper.mle)
}
extendedFamily.absTd <- function(...){extendedFamily.Td(..., family.fun = Family.absTd)}
extendedFamily.Chisq <- function(...){extendedFamily.Td(..., family.fun = Family.Chisq)}
extendedFamily.absTdMixture <- function(...){extendedFamily.TdMixture(..., family.fun = Family.absTdMixture)}
extendedFamily.ChisqMixture <- function(...){extendedFamily.TdMixture(..., family.fun = Family.ChisqMixture)}
#es.extendedFamily.absTd <- function(es, ...)
#{
#	assert.is(es, "xprnSetObject")
#	extendedFamily.absTd(es, df = , ...)
#}
#MLE.absTdMixture <- function(df, ...)
#{
#	stop("MLE.absTdMixture not implemented")
#}
##====================================
Lik.Binom <- function(x, size, prob = 0.5)
{
	dbinom(x = x, size = size, prob = prob, log = FALSE)
}

divergence <- function(P1, P2, base = exp(1), nsample, n, p = 1e-3, ...) # Kullback-Leibler divergence or relative entropy
{
	sca <- if(are(list(P1, P2), "numeric") && missing(nsample) && missing(n))
	{
		stopifnot(length(P1) == length(P2))
		prob.ok <- function(P){is.prob(P = P, ...)}
		stopifnot(prob.ok(P1) && prob.ok(P2))
		lg <- function(x){logb(x = x, base = base)}
		fine <- all(P2 >= 0 | P1 == 0)
		if(!fine)
		{ message("Kullback-Leibler divergence is undefined for those mass functions"); browser()}
		E(x = ifelse(P1 == 0, 0, lg(P1) - lg(P2)), P = P1)
	}
	else if(is(P1, "list") && missing(nsample) && missing(n))
	{
		P1 <- do.call(DensityFunNML, P1)#XXX|:no visible binding for global variable ÔDensityFunNMLÕ
		stopifnot(is(P1, "DensityFunNML"))
		divergence(P1 = P1, P2 = P2, base = base, p = p, ...)
	}
	else if(is(P1, "DensityFunNML") && is(P2, "abscontDistribution") && missing(nsample) && missing(n))
	{
		lims <- integration.limits(object = P1, p = p)
		lower.x <- lims["lower"]
		upper.x <- lims["upper"]
		vector.integrand <- function(X)
		{
			sapply(X, function(x)
			{
				P1(x) * lik.ratio(x = x, object1 = P1, object2 = P2, base = base)
			})
		}
		div <- try(integrate(f = vector.integrand, lower = lower.x, upper = upper.x, ...))
		if(is(div, "try-error"))
		{ message("bad div"); browser()}
		div$value
	}
	else if(are(list(P1, P2), "abscontDistribution") && are(list(n, nsample), "numeric"))
	{
		nsample <- Scalar(nsample)
		rdata1 <- function()
		{
			r(P1)(n = n)
		}
		liks <- sapply(1:nsample, function(i)
		{
			lik.ratio(x = rdata1(), object1 = P1, object2 = P2, base = base)
		})
		assert.is(liks, "numeric")
		stopifnot(length(liks) == nsample)
		mean(liks)
	}
	else
		stop("bad divergence args")
	if(sca < 0)
	{
		warning(paste("divergence approximated as ", sca, " due to insufficient nsample; returning 0"))
		sca <- 0
	}
	ok <- length(sca) == 1 && sca >= 0
	if(!ok)
	{ message("bad sca:"); print(list(sca = sca, P1 = P1, P2 = P2)); browser()}
	Scalar(sca)
}
Divergence <- function(P1, P2, base = exp(1), ...) # see expRelativeEntropy of odd.s, divergence.Bernoulli of data.s, and relative.entropy of verisimilitude.s
{
	base ^ divergence(P1 = P1, P2 = P2, base = base, ...)
}
divergence.Binom <- function(size, min.prob, max.prob, prob, ...)
{
	assert.are(list(size, prob), "numeric")
	stopifnot(length(size) == 1 && size >= 1)
	stopifnot(scalar(min.prob) >= 0 && scalar(max.prob) <= 1 && min.prob <= max.prob)
	stopifnot(length(prob) == 1 && is.prob(prob))
	x <- 0:size
	#XXX|:no visible global function definition for ÔNML.BinomÕ
	divergence(P1 = NML.Binom(x = x, size = size, min.prob = min.prob, max.prob = max.prob), P2 = dbinom(x = x, size = size, prob = prob), ...)
}

new.error <- function(object, ann, estimator)
{
	new("error", Numeric(object), ann = ann, estimator = estimator)
}
#XXX|: to MDL.R
#error <- function(object, P0)
#{
#	assert.is(object, "errorPair")
#	P0 <- Scalar(P0)
#	stopifnot(is.prob(P0))
#	err <- object@error0
#	err@.Data <- P0 * object@error0 + (1 - P0) * object@error1
#	stopifnot(validObject(err))
#	assert.is(err, "error")
#	err
#}

#new.errorPair <- function(error0, error1)
#{
#	new("errorPair", error0 = error0, error1 = error1)
#}
#new.errorPairs <- function(object)
#{
#	new("errorPairs", object)
#}
#-------
risk.estimate <- function(object, root = TRUE, P0, relative)
{
	if(is(object, "error") && missing(P0) && missing(relative))
	{
		av <- mean(object)
		sca <- if(is.function(root))
			root(av)
		else if(root)
			sqrt(av)
		else
			av
		Scalar(sca)
	}
	else if(is(object, "errorPair") && length(P0) == 1 && is.prob(P0))
	{
		if(missing(relative))
			relative <- FALSE
		err <- error(object = object, P0 = P0)
		absolute <- risk.estimate(object = err, root = root)
		if(relative)
		{
			squared <- (1 - P0) ^ 2
			denom <- if(root) sqrt(squared) else squared
			absolute / denom
		}
		else
			absolute
	}
	else
		stop("wrong arguments for risk.estimate")
}

setMethod("annotation", signature(object = "errorPair"), function(object){annotation(object@error0)})

posterior.coverage.error <- function(param.value, posterior.distr, p1, p2, verbose = FALSE, save.time = TRUE, include.Dirac = TRUE, null.param.value = 0)
{
	assert.is(param.value, "numeric")
	assert.is(posterior.distr, "univariateDistribution")
	assert.are(list(p1, p2), "numeric")
	p.ok <- function(p){length(p) == 1 && is.prob(p)}
	stopifnot(p.ok(p1) && p.ok(p2))
	stopifnot(p1 <= p2)
	if(save.time)
	{
		stopifnot(length(param.value) == 1)
#		if(include.Dirac && p1 == 0 && is(posterior.distr, "twoDistrMixture") && is(posterior.distr@distr0, "Dirac") && param.value == q.r(posterior.distr@distr0)(0.5))
		if(include.Dirac && param.value == null.param.value && p2 < P0(posterior.distr))
		{
			if(is(posterior.distr, "twoDistrMixture") && is(posterior.distr@distr0, "Dirac") && null.param.value != q.r(posterior.distr@distr0)(0.5))
				stop("wrong null.param.value?")
			if(verbose && (is(posterior.distr, "SteinChisq") || is(posterior.distr, "SteinChi")))
			{ message("reality check"); browser()}
			0
		}
		else
		{
			prob.le.param.value <- p(posterior.distr)(param.value)
			stopifnot(length(prob.le.param.value) == length(param.value))
			ifelse(prob.le.param.value >= p1 & prob.le.param.value <= p2, 0, 1)
		}
	}
	else
	{
		qfun <- function(p)
		{
			fun <- q.r(posterior.distr)
			if(p == 0)
				-Inf
			else if(p == 1)
				Inf
			else if(!is.na(p))
				fun(p = p)
			else
				stop("bad posterior.coverage.error quantile")
		}
		ci <- c(qfun(p = p1), qfun(p = p2))
		stopifnot(all(!is.na(ci)) && ci[1] <= ci[2])
		ifelse(param.value >= ci[1] & param.value <= ci[2], 0, 1)
	}
}
posterior.squared.error <- function(param.value, posterior.distr, size, verbose = FALSE) # posterior.error.fun = posterior.squared.error
{
	assert.is(param.value, "numeric")
	assert.is(posterior.distr, "univariateDistribution")
	stopifnot(length(size) == 1 && size >= 1 && size == floor(size))
	posterior.param <- r(posterior.distr)(n = size)
	stopifnot(length(param.value) == 1)
	if(verbose)
	{
		print(c(posterior.param = posterior.param, param.value = param.value, err = (posterior.param - param.value) ^ 2))
		browser()
	}
	mean((posterior.param - param.value) ^ 2) # squared.error
}
Distr.error <- function(param.distr, x.distr.fun, x.distr.estimator = NULL, posterior.distr.fun = NULL, posterior.error.fun = NULL, nsample, size = 1, p1 = NULL, p2 = NULL, verbose = FALSE, call.browser = FALSE, ann = date(), ...) # p1 = 1 / 40, p2 = 1 - p1
{
	assert.is(param.distr, "univariateDistribution")
	assert.is(nsample, "numeric")
	if(is.null(p1) && is.null(p2))
	{
		assert.is(size, "numeric")
		stopifnot(length(size) == 1 && size >= 1 && size == floor(size))
	}
	else if(are(list(p1, p2), "numeric"))
	{
		p.ok <- function(p){length(p) == 1 && is.prob(p)}
		stopifnot(p.ok(p1) && p.ok(p2))
		stopifnot(p1 <= p2)
	}
	else
		stop("bad Distr.error args!")
	assert.is(x.distr.fun, "function")
	get.x.distr <- function(param.value)
	{
		assert.is(param.value, "numeric")
		x.distr <- x.distr.fun(param.value)
		assert.is(x.distr, "univariateDistribution")
		x.distr
	}
	get.x <- function(x.distr){assert.is(x.distr, "univariateDistribution"); r(x.distr)(n = size)}
	get.loss <- if(is.null(x.distr.estimator) && !is.null(posterior.distr.fun) && is.function(posterior.error.fun))
	{
		if(is(posterior.distr.fun, "extendedFamily"))
		{
			fam <- posterior.distr.fun
			posterior.distr.fun <- function(x, ...)
			{
				posterior.Distr(object = fam, x = x, ...)
			}
		}
		assert.is(posterior.distr.fun, "function")
		estimator <- posterior.distr.fun
		function(param.value)
		{
			assert.is(param.value, "numeric")
			if(call.browser)
			{ message("param.value for squared error: ", param.value); browser()}
			x.distr <- get.x.distr(param.value)
			x <- get.x(x.distr)
			posterior.distr <- posterior.distr.fun(x = x, call.browser = call.browser, verbose = verbose, ...)
			if(are(list(p1, p2), "numeric"))
				posterior.error.fun(param.value = param.value, posterior.distr = posterior.distr, p1 = p1, p2 = p2, verbose = verbose)
			else if(is.null(p1) && is.null(p2))
				posterior.error.fun(param.value = param.value, posterior.distr = posterior.distr, size = size, verbose = verbose) # posterior.squared.error
			else
				stop("bad posterior.error.fun args!")
		}
	}
	else if(!is.null(x.distr.estimator) && is.null(posterior.distr.fun) && is.null(posterior.error.fun) && is.null(p1) && is.null(p2))
	{
		if(is(x.distr.estimator, "extendedFamily"))
		{
			fam <- x.distr.estimator
			x.distr.estimator <- Distr.estimator(object = fam)
		}
		assert.is(x.distr.estimator, "function")
		estimator <- x.distr.estimator
		function(param.value)
		{
			assert.is(param.value, "numeric")
			if(call.browser)
			{ message("param.value for logarithmic error: ", param.value); browser()}
			x.distr <- get.x.distr(param.value)
			x <- get.x(x.distr)
			x.distr.estimate <- x.distr.estimator(x = x, call.browser = call.browser, verbose = verbose, ...)
			assert.is(x.distr.estimate, "univariateDistribution")
			divergence(P1 = x.distr, P2 = x.distr.estimate, n = 1, nsample = 1) # log_error
		}
	}
	else
		stop("bad args of Distr.error")
	assert.are(list(get.loss, estimator), "function")
	stopifnot(length(nsample) == 1)
	param.values <- r(param.distr)(n = nsample)
	losses <- sapply(param.values, function(param.value)
	{
		assert.is(param.value, "numeric")
		get.loss(param.value = param.value)
	})
	assert.is(losses, "numeric")
	stopifnot(length(losses) == nsample)
	new.error(losses, ann = ann, estimator = estimator)
} # end Distr.error
log_error <- function(param.distr, x.distr.fun, x.distr.estimator = NULL, fam = NULL, nsample, ...)
{
	assert.is(param.distr, "univariateDistribution")
	assert.is(x.distr.fun, "function")
	assert.is(nsample, "numeric")
	if(!is.null(fam))
	{
		assert.is(fam, "extendedFamily")
		stopifnot(is.null(x.distr.estimator))
		message("argument fam is deprecated")
		x.distr.estimator <- fam
	}
	if(is.null(x.distr.estimator))
		stop("x.distr.estimator cannot be created from args")
	Distr.error(param.distr = param.distr, x.distr.fun = x.distr.fun, x.distr.estimator = x.distr.estimator, nsample = nsample, ...)
} # end log_error
squared.error <- function(param.distr, x.distr.fun, posterior.distr.fun, nsample, posterior.error.fun = posterior.squared.error, ...)
{# posterior.error.fun = posterior.coverage.error of coverage.errorPairs.absTd will not give squared error

	assert.is(param.distr, "univariateDistribution")
	assert.is(x.distr.fun, "function")
	assert.is(posterior.distr.fun, c("function", "extendedFamily"))
	assert.is(nsample, "numeric")
	Distr.error(param.distr = param.distr, x.distr.fun = x.distr.fun, posterior.distr.fun = posterior.distr.fun, posterior.error.fun = posterior.error.fun, nsample = nsample, ...)
} # end squared.error
squared.error.noncentralDistr <- function(P0, df, ncp, posterior.distr.fun = NULL, lower.P0, lower.ncp, ann = NULL, Distr.fun = stop("try something like absTd"), family.fun = stop("try something like extendedFamily.absTdMixture"), ...)
{
	assert.is(Distr.fun, "function")
	stopifnot(is.prob(P0) && length(P0) == 1)
	df <- Scalar(df)
	ncp <- scalar(ncp)
	if(is.null(posterior.distr.fun))
	{
		ann0 <- "confidence posterior"
		posterior.distr.fun <- if(missing(lower.P0) && missing(lower.ncp))
		{
			stop("generic confidence posterior not supported")
		}
		else if(TRUE)
		{
			if(missing(lower.P0))
				lower.P0 <- 0
			if(missing(lower.ncp))
				lower.ncp <- 1
			ann0 <- paste("LFDR; [", lower.P0, ",1]",  sep = "")
			family.fun(df = df, lower.P0 = lower.P0, lower.ncp = lower.ncp)
		}
		else
			stop("bad arguments for squared.error.noncentralDistr")
	}
	else
		ann0 <- date()
	if(is.null(ann))
		ann <- ann0
	param.distr <- new.twoDistrMixture(P0 = P0, distr0 = Dirac(location = 0), distr1 = Dirac(location = ncp))
	x.distr.fun <- function(NCP){Distr.fun(df = df, ncp = NCP)}
	squared.error(param.distr = param.distr, x.distr.fun = x.distr.fun, posterior.distr.fun = posterior.distr.fun, ann = ann, ...) # P0 = P0, distr0 = get.distr(NCP = 0), distr1 = get.distr(NCP = ncp)
} # end squared.error.noncentralDistr
squared.error.absTd <- function(df = Inf, posterior.distr.fun = NULL, lower.P0, lower.ncp, ann = NULL, shrink = TRUE, ...)
{
	assert.is(df, "numeric")
	if(is.null(posterior.distr.fun))
	{
		if(missing(lower.P0) && missing(lower.ncp))
		{
			posterior.distr.fun <- if(df == Inf)
			{
				function(x, ...)
				{
					Stein.posterior.Distr(x = x, shrink = shrink, squared.parameter = FALSE, ...)
				}
			}
			else
				stop("not yet implemented: confidence posterior with finite df")
			if(is.null(ann))
			{
				ann <- if(shrink)
					"marginal confidence" # shrinks toward null distribution
				else
					"conditional confidence" # conditional on alternative distribution
			}
		} # end confidence posterior
	}
	else
		stopifnot(shrink)
	squared.error.noncentralDistr(df = df, posterior.distr.fun = posterior.distr.fun, lower.P0 = lower.P0, lower.ncp = lower.ncp, Distr.fun = absTd, family.fun = extendedFamily.absTdMixture, ann = ann, ...)
} # end squared.error.absTd
squared.error.Chisq <- function(df = 1, posterior.distr.fun = NULL, lower.P0, lower.ncp, ann = NULL, ...)
{
	assert.is(df, "numeric")
	if(is.null(posterior.distr.fun))
	{
		if(missing(lower.P0) && missing(lower.ncp))
		{
			posterior.distr.fun <- if(df == 1)
			{
				function(x, ...)
				{
					SteinChisq(ncp = x, df = df) # Stein.posterior.Distr(x = x, shrink = shrink, squared.parameter = TRUE, ...) has wrong x
				}
			}
			else
				stop("not yet tested: confidence posterior with non-unit df")
			if(is.null(ann))
			{
				ann <- "confidence" # shrinks toward null distribution
			}
		} # end confidence posterior
	}
	squared.error.noncentralDistr(df = df, posterior.distr.fun = posterior.distr.fun, lower.P0 = lower.P0, lower.ncp = lower.ncp, Distr.fun = Chisq, family.fun = extendedFamily.ChisqMixture, ann = ann, ...)
} # end squared.error.Chisq
squared.errorPair <- function(FUN = squared.error.noncentralDistr, ...)
{
	get.error <- function(P0){FUN(P0 = P0, ...)}
	error0 <- get.error(P0 = 1) # null
	error1 <- get.error(P0 = 0) # alternative
	new.errorPair(error0 = error0, error1 = error1)
}
squared.errorPair.absTd <- function(...)
{
	squared.errorPair(FUN = squared.error.absTd, ...)
}
squared.errorPair.Chisq <- function(...)
{
	squared.errorPair(FUN = squared.error.Chisq, ...)
}
squared.errorPairs <- function(arg.lists, nsample = stop("missing nsample"), seed = numeric(0), FUN = squared.errorPair, ...)
{
	assert.is(arg.lists, "list")
	stopifnot(length(arg.lists) >= 1)
	lis <- lapply(arg.lists, function(arg.list)
	{
		assert.is(arg.list, "list")
		setSeed(seed = seed)
		pair <- do.call(FUN, c(arg.list, list(..., nsample = nsample)))
		message("  completed ", annotation(pair), " on ", date(), sep = "")
		pair
	})
	names(lis) <- names(arg.lists)
	new.errorPairs(lis)
}
squared.errorPairs.absTd <- function(arg.lists = stop("missing arg.lists"), ncp = stop("missing ncp"), ...)
{
	squared.errorPairs(arg.lists = arg.lists, FUN = squared.errorPair.absTd, ncp = ncp, ...)
}
coverage.errorPairs.absTd <- function(arg.lists = stop("missing arg.lists"), ncp = stop("missing ncp"), p1 = 0, p2, ...)
{
	assert.are(list(p1, p2), "numeric")
	squared.errorPairs.absTd(arg.lists = arg.lists, ncp = ncp, posterior.error.fun = posterior.coverage.error, p1 = p1, p2 = p2, ...)
}
squared.errorPairs.Chisq <- function(arg.lists = stop("missing arg.lists"), ncp = stop("missing ncp"), ...)
{
	squared.errorPairs(arg.lists = arg.lists, FUN = squared.errorPair.Chisq, ncp = ncp, ...)
}

quadratic.error <- function(indicator0, distr0, distr1, nsample = stop("nsample missing"), indicator0.estimator = NULL, size = 1, verbose = FALSE, call.browser = FALSE, alpha, ann = date(), ...)
{
	if(is.logical(indicator0))
		indicator0 <- as.numeric(indicator0)
	assert.is(indicator0, "numeric")
	stopifnot(all(indicator0 %in% c(0, 1)))
	assert.is(nsample, "numeric")
	if(missing(nsample) && length(indicator0) > 1)
		nsample <- length(indicator0)
	stopifnot(length(nsample) == 1)
	if(length(indicator0) == 1)
		indicator0 <- rep(indicator0, nsample)
	stopifnot(nsample == length(indicator0))
	if(is(indicator0.estimator, "extendedFamily"))
	{
		fam <- indicator0.estimator
		indicator0.estimator <- function(x)
		{
			# mle.out <- mle(object = fam, x = x, call.browser = call.browser, verbose = verbose, ...)
			posteriorP0(distr0 = distr0, fam = fam, x = x, call.browser = call.browser, verbose = verbose, ...) # P0(mle.out)
		}
	}
	assert.is(indicator0.estimator, "function")
	if(!missing(alpha) && !is.nothing(alpha))
	{
		stopifnot(length(alpha) == 1 && is.prob(alpha))
		ann <- paste(ann, "; level: ", alpha, sep = "")
		P0.fun <- indicator0.estimator
		indicator0.estimator <- function(x)
		{
			p0 <- P0.fun(x)
			if(p0 <= alpha) 0 else 1
		}
	}
	losses <- sapply(indicator0, function(null.indicator)
	{
		if(call.browser)
		{ message("null.indicator: ", null.indicator); browser()}
		x.distr <- if(null.indicator == 1)
			distr0
		else if(null.indicator == 0)
			distr1
		else
			stop("bad null.indicator")
		assert.is(x.distr, "univariateDistribution")
		x <- r(x.distr)(n = size)
		estimate0 <- indicator0.estimator(x)
		assert.is(estimate0, "numeric")
		stopifnot(length(estimate0) == 1 && is.prob(estimate0))
		(estimate0 - null.indicator) ^ 2
	})
	assert.is(losses, "numeric")
	stopifnot(length(losses) == nsample)
	new.error(losses, ann = ann, estimator = indicator0.estimator)
} # end quadratic.error
quadratic.error.noncentralDistr <- function(indicator0, df, ncp, indicator0.estimator = NULL, lower.P0, lower.ncp, ann, Distr.fun = stop("try something like absTd"), family.fun = stop("try something like extendedFamily.absTdMixture"), ...)
{
	df <- Scalar(df)
	ncp <- scalar(ncp)
	if(is.null(indicator0.estimator))
	{
		ann0 <- "p-value"
		indicator0.estimator <- if(missing(lower.P0) && missing(lower.ncp))
		{
			function(x)
			{
				posteriorP0(distr0 = Distr.fun(df = df, ncp = 0), x = x) # changed 110318 without testing.
			}
		}
		else if(TRUE) # missing(alpha))
		{
			if(missing(lower.P0))
				lower.P0 <- 0
			if(missing(lower.ncp))
				lower.ncp <- 1
			ann0 <- paste("LFDR; [", lower.P0, ",1]",  sep = "")
			family.fun(df = df, lower.P0 = lower.P0, lower.ncp = lower.ncp)
		}
		else
			stop("bad arguments for quadratic.error.noncentralDistr")
	}
	else
		ann0 <- date()
	if(missing(ann))
		ann <- ann0
	get.distr <- function(NCP){Distr.fun(df = df, ncp = NCP)}
	quadratic.error(indicator0 = indicator0, distr0 = get.distr(NCP = 0), distr1 = get.distr(NCP = ncp), indicator0.estimator = indicator0.estimator, ann = ann, ...)
}
quadratic.error.absTd <- function(...)
{
	quadratic.error.noncentralDistr(Distr.fun = absTd, family.fun = extendedFamily.absTdMixture, ...)
}
quadratic.error.Chisq <- function(...)
{
	quadratic.error.noncentralDistr(Distr.fun = Chisq, family.fun = extendedFamily.ChisqMixture, ...)
}
quadratic.errorPair <- function(FUN = quadratic.error, ...)
{
	get.error <- function(indicator0){FUN(indicator0 = indicator0, ...)}
	error0 <- get.error(indicator0 = 1) # null
	error1 <- get.error(indicator0 = 0) # alternative
	new.errorPair(error0 = error0, error1 = error1)
}
quadratic.errorPair.absTd <- function(...)
{
	quadratic.errorPair(FUN = quadratic.error.absTd, ...)
}
quadratic.errorPair.Chisq <- function(...)
{
	quadratic.errorPair(FUN = quadratic.error.Chisq, ...)
}
quadratic.errorPairs <- function(arg.lists, nsample = stop("missing nsample"), seed = numeric(0), FUN = quadratic.errorPair, ...)
{
	assert.is(arg.lists, "list")
	stopifnot(length(arg.lists) >= 1)
	lis <- lapply(arg.lists, function(arg.list)
	{
		assert.is(arg.list, "list")
		setSeed(seed = seed)
		pair <- do.call(FUN, c(arg.list, list(..., nsample = nsample)))
		message("  completed ", annotation(pair), " on ", date(), sep = "")
		pair
	})
	names(lis) <- names(arg.lists)
	new.errorPairs(lis)
}
quadratic.errorPairs.absTd <- function(arg.lists = stop("missing arg.lists"), df = stop("missing df"), ncp = stop("missing ncp"), ...)
{
	quadratic.errorPairs(arg.lists = arg.lists, FUN = quadratic.errorPair.absTd, df = df, ncp = ncp, ...)
}
quadratic.errorPairs.Chisq <- function(arg.lists = stop("missing arg.lists"), df = stop("missing df"), ncp = stop("missing ncp"), ...)
{
	quadratic.errorPairs(arg.lists = arg.lists, FUN = quadratic.errorPair.Chisq, df = df, ncp = ncp, ...)
}

Stein.posterior.Distr <- function(x, shrink = TRUE, call.browser = FALSE, verbose = FALSE, df = Inf, squared.parameter = FALSE)
{
	if(call.browser)
	{ message("you called?"); browser()}
	assert.is(x, "numeric")
	stopifnot(all(x >= 0))
	stopifnot(length(df) == 1 && df == Inf)
	dis <- if(shrink)
	{
		get.distr <- if(squared.parameter) SteinChisq else SteinChi
		get.distr(ncp = sum(x ^ 2), df = length(x)) # confidence distribution from distribution of abs(x)
	}
	else if(length(x) == 1 && !squared.parameter)
		absTd(ncp = x, df = df) # no shrinkage; ncp^2 has chi^2 distr with ncp=x
	else
		stop("bad Stein.posterior.Distr args")
	if(verbose)
		print(dis)
	stopifnot(is(dis, "univariateDistribution"))
	dis
}
posterior.Distr <- function(object, x, unknown.param.name, null.param.value = 0, ...)
{
	assert.is(x, "numeric")
	dis <- if(is(object, "extendedFamily"))
	{
		mle.out <- mle(object = object, x = x, ...)
		assert.is(mle.out, "nmle")#XXX|:changing from mle to nmle
		prior0 <- P0(mle.out)
		ok <- is.prob(prior0) && length(prior0) == 1
		if(!ok)
		{ message("bad prior0 in posterior.Distr: ", prior0); browser()}
		unknown.param <- mle.out@unknown.param
		stopifnot(length(unknown.param) >= 1)
		wrong.param.name <- "P0"
		if(missing(unknown.param.name))
		{
			nam <- setdiff(names(unknown.param), wrong.param.name)
			unknown.param.name <- if(length(nam) == 1)
				nam
			else
			{
				message("parameter of interest was not specified")
				browser()
			}
		}
		stopifnot(length(unknown.param.name) == 1)
		if(unknown.param.name == wrong.param.name)
			stop("incorrect parameter of interest; see posteriorP0")
		unknown.param <- unknown.param[unknown.param.name]
		stopifnot(length(unknown.param) == 1)
		alt.param.value <- unknown.param[[1]]
		assert.is(alt.param.value, "numeric")
		stopifnot(length(alt.param.value) == 1)
		new.twoDistrMixture(P0 = prior0, distr0 = Dirac(location = null.param.value), distr1 = Dirac(location = alt.param.value))
	}
	else
		stop("posterior.Distr not yet implemented for object of that class")
	assert.is(dis, "univariateDistribution")
	dis
}
Distr.estimate <- function(object, x, ...)
{
	assert.is(object, "extendedFamily")
	assert.is(x, "numeric")
	mle.out <- mle(object = object, x = x, ...)
	assert.is(mle.out, "nmle")#XXX|:changing from mle to nmle
	prior0 <- P0(mle.out)
	ok <- is.prob(prior0) && length(prior0) == 1
	if(!ok)
	{ message("bad prior0 in Distr.estimate: ", prior0); browser()}
	distr1.estimate <- Distr(object = object, unknown.param = mle.out@unknown.param)
	assert.is(distr1.estimate, "univariateDistribution")
	if(FALSE) # potential problem found 7 March 2011, but ok if object is mixture (prior0 is for error checking)
	{
		message("Distr.estimate will return distr1.estimate, ignoring prior0 of ", prior0)
		browser()
	}
	distr1.estimate
}
Distr.estimator <- function(object)
{
	function(x, ...)
	{
		Distr.estimate(object = object, x = x, ...)
	}
}

new.posteriorP0 <- function(object, estimate = new.mle(), x = numeric(0), distr0 = NULL, family = NULL)
{
	if(is.null(distr0))
		distr0 <- blank.Distr()
	if(is.null(family))
		family <- blank.Family()
	post0 <- Numeric(object)
	if(is.null(names(x)) && length(x) == length(post0))
		names(x) <- names(post0)
	if(!sameNames(post0, x))
	{ message("name problem"); browser()}
	new("posteriorP0", post0, estimate = estimate, x = x, distr0 = distr0, family = family)
}
posteriorP0 <- function(distr0 = NULL, fam = NULL, x, verbose = FALSE, individual, confidence, unknown.param = NULL, save.time = logical(0), ...) # save.time = logical(0) on 110417; was FALSE
{# confidence in null hypothesis or estimated local false discovery rate

	if(is(distr0, "blank.Distr"))
		distr0 <- NULL
	pP0 <- if(!missing(individual) && individual && is(x, "numeric") && length(x) > 1)
	{
		vec <- sapply(x, function(X){posteriorP0(distr0 = distr0, fam = fam, x = X, verbose = verbose, confidence = confidence, unknown.param = unknown.param, save.time = save.time, ...)})
		assert.is(vec, "numeric")
		stopifnot(length(vec) == length(x))
		names(vec) <- names(x)
		stopifnot(is.prob(vec))
		new.posteriorP0(vec, x = x, family = fam, distr0 = distr0)
	}
	else if(is(distr0, "univariateDistribution") && is.null(fam) && is(x, "numeric") && length(x) == 1 && is.nothing(list(...)))
	{
		if(missing(confidence))
			confidence <- TRUE
		stopifnot(confidence)
		pval <- p(distr0)(q = x, lower.tail = FALSE)
		stopifnot(length(pval) == 1 && is.prob(pval))
		new.posteriorP0(pval, x = x, family = fam, distr0 = distr0)
	}
	else if(is.null(fam) && is(x, "xprnSetObject"))
	{
		df <- Df(x)
		if(is.null(distr0))
			distr0 <- absTd(df = df, ncp = 0)
		assert.is(distr0, "univariateDistribution")
		if(missing(confidence))
			confidence <- is.nothing(list(...))
		if(missing(individual))
			individual <- confidence
		at <- stat(x, FUN = statistic.absTd)
		if(is.null(fam) && !confidence)
			fam <- extendedFamily.absTdMixture(df = df, save.time = save.time, ...) # lower.ncp = 1 / 1e3, lower.P0 = 90 / 100
		if(individual)
		{
			vec <- posteriorP0(distr0 = distr0, fam = fam, x = at, verbose = verbose, individual = individual, confidence = confidence, unknown.param = unknown.param, save.time = save.time)
			stopifnot(length(vec) == Length(x))
			names(vec) <- names(x)
			stopifnot(is.prob(vec))
			vec
		}
		else if(!confidence) # default
		{
			assert.is(fam, "extendedFamily")
			posteriorP0(distr0 = distr0, fam = fam, x = at, verbose = verbose, individual = individual, confidence = confidence, unknown.param = unknown.param, save.time = save.time)
		}
		else
			stop("posteriorP0 with xprnSetObject has not yet been implemented for those arguments")
	}
	else if(is(fam, "Family") && is(x, "numeric") && length(x) >= 1)
	{
		if(missing(individual))
			individual <- FALSE
		if(individual)
			stop("individual should have been caught earlier")
		if(missing(confidence))
			confidence <- FALSE
		if(confidence)
			stop("confidence should have been caught earlier")
		prior0 <- if(is.null(unknown.param))
		{
			assert.is(fam, "extendedFamily")
			mle.out <- mle(object = fam, x = x, ...)
			assert.is(mle.out, "nmle")#XXX|:changing from mle to nmle
			unknown.param <- mle.out@unknown.param
			P0(mle.out)
		}
		else
		{
			unknown.param <- unknownParam(unknown.param)
			mle.out <- new.mle(unknown.param = unknown.param)
			unknown.param$P0
		}
		if(prior0 > 1)
		{
			message("    prior0 of ", round(prior0, 2), " reset to 1 on ", date())
			prior0 <- 1
		}
		if(prior0 < 0)
		{
			message("    prior0 of ", round(prior0, 2), " reset to 0 on ", date())
			prior0 <- 0
		}
		ok <- is.prob(prior0) && length(prior0) == 1
		if(!ok)
		{ message("bad prior0 in posteriorP0: ", prior0); browser()}
		if("P0" %in% names(unknown.param))
		{
			unknown.param$P0 <- prior0
			mle.out@unknown.param <- unknown.param
		}
		else
			stop("cannot correct unknown.param")
		distr.estimate <- Distr(object = fam, unknown.param = unknown.param) # mixture by default; error corrected 110413
		assert.are(list(distr0, distr.estimate), "univariateDistribution")
		get.dens <- function(distr){sapply(x, d(distr))} # (x)
		dens0 <- get.dens(distr0)
		dens <- get.dens(distr.estimate)
		if(verbose)
		{
			print(unknown.param)
			print(c(prior0 = prior0, dens0 = dens0, dens = dens))
		}
		post0 <- prior0 * dens0 / dens # (prior0 * dens0 + (1 - prior0) * dens1)
		stopifnot(all(length(post0) == c(length(dens0), length(dens), length(x))))
		nan.boo <- is.nan(post0)
		if(any(nan.boo))
		{
			replacement <- if(all(nan.boo)) prior0 else median(post0[!nan.boo])
			stopifnot(length(replacement) == 1)
			message("\n    ", sum(nan.boo), " of ", length(nan.boo), " LFDR estimates replaced with ", replacement)
			post0 <- ifelse(nan.boo, replacement, prior0)
			print(stats(post0))
			cat("\n\n")
		}
		if(!is.prob(post0))
		{
			message("is.prob(post0) is not TRUE; post0: ", post0); browser()
		}
		names(post0) <- names(x)
		if(!sameNames(post0, x))
		{ message("impossible name problem"); browser()}
		try(new.posteriorP0(object = post0, estimate = mle.out, x = x, family = fam, distr0 = distr0))
	}
	else if(missing(individual) && missing(confidence))
		stop("posteriorP0 has not yet been implemented for those arguments")
	else
	{	message("posteriorP0 has not yet been implemented for those arguments"); browser()} # end pP0
	if(is.err(pP0))
	{ message("name problem?"); browser()}
	assert.is(pP0, "posteriorP0")
	pP0
} # end posteriorP0

parametric.bootstrap.estimate <- function(object, estimator, fam, unknown.param, n = NULL)
{
	if(!missing(object) && missing(fam) && missing(unknown.param))
	{
		if(is(object, "resampled.posteriorP0") && missing(estimator))
		{
			estimator <- object@posteriorP0.fun
			fam <- object@original@family
			unknown.param <- object
			if(is.null(n))
				n <- length(object@original)
		}
		else if(is(object, "posteriorP0") && is.function(estimator))
		{
			fam <- object@family
			unknown.param <- object
			if(is.null(n))
				n <- length(object)
		}
		else
			stop("bad arguments in parametric.bootstrap.estimate........")
		parametric.bootstrap.estimate(estimator = estimator, fam = fam, unknown.param = unknown.param, n = n)
	}
	else if(missing(object) && !missing(estimator) && !missing(fam) && !missing(unknown.param) && !missing(n))
	{
		assert.is(estimator, "function")
		assert.is(fam, "Family")
		assert.is(n, "numeric")
		stopifnot(length(n) == 1 && n >= 1)
		distr <- Distr(object = fam, unknown.param = unknown.param)
		n1 <- floor(n)
		n2 <- ceiling(n)
		if(n1 != n2)
		{
			prob.n2 <- n2 - n1
			stopifnot(is.prob(prob.n2))
			n <- if(runif(n = 1) < prob.n2) n2 else max(n1, 1)
		}
		estimator(r(distr)(n = n))
	}
	else
		stop("bad arguments in parametric.bootstrap.estimate!!!!!!!!")
}
##------------------------
#see http://stat.stanford.edu/~omkar/monograph/simz.R 
myzsim <- function(alpha=0,N=1000,n=1,stand=1,J=5){c(alpha,N,n,stand,J)}
#------------------------
r.xprnSet <- function(alpha = default(0, "alpha"), nfeature, size, mean1 = NULL, stand = default(0, "stand", verbose = FALSE), J)
{
	stopifnot(is.prob(alpha))
	assert.is(nfeature, "numeric")
	assert.is(size, "numeric")
	N <- nfeature
	n <- size
	if(missing(J))
	{
		J <- 5
		ratio <- floor(N / J)
		is.ok <- function(J){ratio >= 1 && ratio == N / J}
		if(!is.ok(J))
			J <- N
		stopifnot(is.ok(J))
		default(J, "J", verbose = FALSE)
	}
	exprs <- myzsim(alpha = alpha, N = N, n = n, stand = stand, J = J)
	#XXX|:no visible global function definition for ÔmyzsimÕ
	assert.is(exprs, "matrix")
	if(!is.null(mean1))
	{
		assert.is(mean1, "numeric")
		stopifnot(length(mean1) <= nrow(exprs))
		for(i in 1:length(mean1))
			exprs[i, ] <- exprs[i, ] + mean1[i]
	}
	stopifnot(size == ncol(exprs))
	es <- xprnSet(exprs = exprs)
	stopifnot(nfeature == Length(es))
	es
}
r.xprnSetPair <- function(nfeature, x.size, y.size = NULL, x.mean1 = NULL, y.mean1 = NULL, ...)
{
	if(is.null(y.size))
		y.size <- x.size
	assert.are(list(nfeature, x.size, y.size), "numeric")
	get.xprnSet <- function(size, mean1)
	{
		r.xprnSet(nfeature = nfeature, size = size, mean1 = mean1, ...)
	}
	xprnSetPair(x = get.xprnSet(size = x.size, mean1 = x.mean1), y = get.xprnSet(size = y.size, mean1 = y.mean1))
}
r.xprnSetObject <- function(alpha = default(0, "alpha"), nfeature, size = NULL, mean1 = NULL, x.size = NULL, y.size = NULL, x.mean1 = NULL, y.mean1 = NULL, J = NULL, ...)
{
	assert.is(alpha, "numeric")
	assert.is(nfeature, "numeric")
	if(is.null(J))
		J <- 5
	fun <- if(is.null(size) && is.null(x.size))
		stop("size and x.size are both missing")
	else if(are.null(list(size, mean1)))
	{
		assert.is(x.size, "numeric")
		function(...){r.xprnSetPair(..., x.size = x.size, y.size = y.size, x.mean1 = x.mean1, y.mean1 = y.mean1)}
	}
	else if(are.null(list(x.size, y.size, x.mean1, y.mean1)))
	{
		assert.is(size, "numeric")
		function(...){r.xprnSet(..., size = size, mean1 = mean1)}
	}
	else
		stop("bad size or x.size")
	fun(alpha = alpha, nfeature = nfeature, J = J, ...)
}
#-----------------
# this alpha program (Efron, 2010, pp. 148, 249-250) estimates the root mean square correlation by Efron's equation (8.18)
myalpha <- function( dat , samp=10,start=1, n = ncol(dat),perm=0, sw = 0){c(dat , samp,start, n = ncol(dat),perm, sw)}

##---------------
RMS.cor <- function(object, samp = -1, start = 1, n = NULL, location.fun = mean, scale.fun = NULL)
{
	recall <- function(dat)
	{
		RMS.cor(object = dat, samp = samp, start = start, n = n, location.fun = location.fun, scale.fun = scale.fun)
	}
	if(is(object, "XprnSetObject") || is(object, "XprnSetPair"))
		recall(logb(object))
	else if(is(object, "xprnSetPair"))
		sqrt(mean(c(recall(object@x), recall(object@y)) ^ 2))
	else if(is(object, "xprnSet"))
		recall(exprs(object))
	else if(is(object, "matrix"))
	{
		dat <- standardize(object, location.fun = location.fun, scale.fun = scale.fun) # added 110412
		a <- myalpha(dat = dat, samp = samp, start = start, n = if(is.null(n)) ncol(object) else n) # in Efron.s
		#XXX|no visible global function definition for ÔmyalphaÕ
		stopifnot(is.prob(a))
		Scalar(a)
	}
	else
		stop("incorrect arguments for RMS.cor")
}
effective.Length <- function(object, verbose = FALSE, ...) # Efron, AOAS 2009, (4.8)
{
	assert.is(object, "xprnSetObject")
	a <- RMS.cor(object = object, ...)
	if(verbose)
		message("alpha from RMS.cor: ", round(a, 3))
	N <- Length(object)
	divisor <- 1 + (N - 1) * a ^ 2
	stopifnot(divisor <= N && divisor >= 1)
	Scalar(N / divisor)
}

new.resampled.posteriorP0 <- function(object, unknown.params, x, original, posteriorP0.fun)
{
	assert.is(object, "matrix")
	new("resampled.posteriorP0", object, unknown.params = new.unknownParams(unknown.params), x = x, original = original, posteriorP0.fun = posteriorP0.fun)
}
resampled.posteriorP0 <- function(object = posteriorP0, x, nsample, n = NULL, individual = FALSE, confidence = FALSE, lower.ncp = NULL, lower.P0 = NULL, upper.P0 = NULL, assumed.P0 = NULL, verbose = TRUE, ...)
{
	assert.is(nsample, "numeric")
	assert.is(object, "function") # returns posteriorP0 object
	if(is.null(n) && is(x, "xprnSetObject")) # added 110410
	{
		n <- default(effective.Length(object = x, verbose = verbose), name = "  n", verbose = verbose)
		if(verbose)
			message(round(n, 2), " effective features; ", Length(x), " actual features")
	}
	get.original <- function(x, ...)
	{
		orig <- try(object(x = x, individual = individual, confidence = confidence, ...))
		if(is.err(orig))
		{ message("bad orig"); print(list(...)); browser()}
		orig
	}
	original <- get.original(x = x, lower.ncp = lower.ncp, lower.P0 = lower.P0, upper.P0 = upper.P0, assumed.P0 = assumed.P0, ...)
	assert.is(original, "posteriorP0")
	posteriorP0.fun <- if(is(x, "numeric"))
	{
		stop("lower.ncp, lower.P0, etc. are not yet handled for this argument")
		get.original
	}
	else if(is(original@family, "extendedFamily") && is(original@distr0, "univariateDistribution") && !is(original@distr0, "blank.Distr"))
	{
		if(confidence || individual)
			stop("the slots are not yet supported")
		function(x, unknown.param = NULL){get.original(x = x, fam = original@family, distr0 = original@distr0, unknown.param = unknown.param, ...)} # lower.ncp and lower.P0 are in original@family
	}
	else
	{ message("cannot get posteriorP0.fun"); browser()}
	nro <- length(original)
	nco <- nsample
	get.mat <- function()
	{
		mat <- matrix(nrow = nro, as.numeric(rep(NA, nro * nco)))
		rownames(mat) <- names(original)
		stopifnot(nco == ncol(mat))
		mat
	}
	posteriorP0.mat <- x.mat <- get.mat()
	unknown.params <- list()
	is.vec.ok <- function(vec){length(vec) == nrow(posteriorP0.mat) && all(names(vec) == rownames(posteriorP0.mat))}
	for(j in 1:nsample)
	{
		pbe <- parametric.bootstrap.estimate(object = original, estimator = posteriorP0.fun, n = n)
		unknown.param <- unknownParam(pbe)
		unknown.params <- c(unknown.params, list(unknown.param))
		posteriorP0.vec <- try(posteriorP0.fun(x = original@x, unknown.param = unknown.param))
		if(is.err(posteriorP0.vec) || !is.vec.ok(posteriorP0.vec))
		{ message("bad posteriorP0.vec"); browser()}
		posteriorP0.mat[, j] <- posteriorP0.vec
		x.vec <- posteriorP0.vec@x
		stopifnot(is.vec.ok(x.vec))
		x.mat[, j] <- x.vec
	}
	stopifnot(all(is.finite(posteriorP0.mat)))
	stopifnot(all(is.finite(x.mat)))
	new.resampled.posteriorP0(object = posteriorP0.mat, unknown.params = unknown.params, x = x.mat, original = original, posteriorP0.fun = posteriorP0.fun)
} # end resampled.posteriorP0
#XXX|: should change all "mean." by "mean_"?
mean.posteriorP0 <- function(..., method = default("p", "mean.posteriorP0 method"))
{
	rpost0 <- resampled.posteriorP0(...)
	mean(rpost0, method = method) # as(rpost0, "posteriorP0")
}
#XXX|: should all "Median." by "Median_"?
Median.posteriorP0 <- function(..., method = default("p", "Median.posteriorP0 method"), nresample = default(20, "Median.posteriorP0 nresample"))
{
	rpost0 <- resampled.posteriorP0(..., nresample = nresample)
	Median(rpost0, method = method) # as(rpost0, "posteriorP0")
}

new.posteriorP0.error <- function(object, portion.above, bias, P0.bias, ncp.bias, score, null.indicator, unknown.param, ann, estimator, x.fun)
{
	err <- new.error(object = object, ann = ann, estimator = estimator)
	new("posteriorP0.error", err, portion.above = Numeric(portion.above), bias = bias, P0.bias = scalar(P0.bias), ncp.bias = scalar(ncp.bias), score = Numeric(score), null.indicator = Numeric(null.indicator), unknown.param = unknown.param, x.fun = x.fun)
}
posteriorP0.error <- function(posteriorP0.fun, true.posteriorP0.fun = posteriorP0, method, nresample = NULL, nsim, P0, alpha = NULL, nfeature = NULL, size = NULL, mean1 = NULL, x.size = NULL, y.size = x.size, x.mean1 = NULL, y.mean1 = NULL, J = NULL, unknown.param, null.indicator, x.fun, ann, min.P0.discrepancy = 0, max.P0.discrepancy = 1, save.time = logical(0), debug.null.indicator = FALSE, ...) #  save.time = logical(0) on 110417; was default(TRUE, "save.time") (save.time added 110415)
{
	if(!missing(P0))
	{
		stopifnot(is.prob(P0) && length(P0) == 1)
		assert.is(nfeature, "numeric")
		get.alt.mean <- function(mean1.arg)
		{
			if(is.null(mean1.arg))
				mean1.arg
			else if(length(mean1.arg) == 1)
			{
				frac <- (1 - P0) * nfeature
				nfeature1 <- round(frac)
				if(abs(frac - nfeature1) > 1e-2)
				{ message("nfeature1: ", nfeature1); browser()}
				rep(mean1.arg, nfeature1)
			}
			else
				stop("bad mean1.arg for given P0")
		}
		if(are.null(list(mean1, x.mean1, y.mean1)))
			stop("P0 cannot be used without specifying one of these: mean1, x.mean1, y.mean1")
		mean1 <- get.alt.mean(mean1)
		x.mean1 <- get.alt.mean(x.mean1)
		y.mean1 <- get.alt.mean(y.mean1)
	}
	if(missing(posteriorP0.fun))
	{
		if(missing(method))
			method <- default("p", "posteriorP0.error method")
		if(is.null(nresample))
		{
			nresample <- if(method %in% c("p", "bc")) # p = percentile; bc = bias-corrected
				1
			else if(method == "bcp") # bcp = bias-corrected percentile
				20
			else
				stop("method unrecognized")
			default(nresample, paste("posteriorP0.error", method, "nresample"))
		}
		posteriorP0.fun <- function(...){mean.posteriorP0(..., method = method, nsample = nresample)}
		if(missing(ann))
			ann <- "mean.posteriorP0"
	}
	else if(missing(method) && is.null(nresample))
	{
		if(missing(ann))
			ann <- "posteriorP0"
	}
	else
		stop("awful arguments")
	assert.is(posteriorP0.fun, "function")
	if(missing(unknown.param) || missing(null.indicator))
	{
		assert.are(list(alpha, nfeature), "numeric")
		get.mean <- function(mean1)
		{
			vec <- c(if(is.null(mean1)) numeric(0) else mean1, rep(0, nfeature - length(mean1)))
			assert.is(vec, "numeric")
			stopifnot(length(vec) == nfeature)
			vec
		}
		diff.means <- if(is.null(mean1) && !missing(x.mean1))
		{
			tsize <- t_size(x.size, y.size)
			get.mean(x.mean1) - get.mean(y.mean1)
		}
		else if(!is.null(mean1) && missing(x.mean1))
		{
			tsize <- size
			get.mean(mean1)
		}
		else
			stop("cannot get unknown.param")
		diff.mean1 <- diff.means[1]
		if(any(diff.means[diff.means != 0] != diff.mean1))
			stop("cannot get ncp for missing unknown.param")
		if(missing(null.indicator))
			null.indicator <- ifelse(diff.means == 0, 1, 0)
		stopifnot(is.character(null.indicator) || length(null.indicator) == nfeature)
		if(missing(unknown.param))
		{
			assert.is(null.indicator, "numeric")
			pi0 <- sum(null.indicator == 1) / nfeature
			if(missing(P0))
				P0 <- pi0
			else
				stopifnot(P0 == pi0)
			unknown.param <- list(P0 = P0, ncp = diff.mean1 * sqrt(tsize))
		}
	}
	unknown.param <- unknownParam(unknown.param)
	assert.is(unknown.param, "list")
	assert.is(nsim, "numeric")
	if(missing(x.fun))
	{
		assert.are(list(alpha, nfeature), "numeric")
		x.fun <- function()
		{
			x <- r.xprnSetObject(alpha = alpha, nfeature = nfeature, size = size, mean1 = mean1, x.size = x.size, y.size = y.size, x.mean1 = x.mean1, y.mean1 = y.mean1, J = J)
			assert.is(x, "xprnSetObject")
			x
		}
	}
	else
		stopifnot(are.null(list(alpha, size, mean1, x.size, y.size, x.mean1, y.mean1, J))) # nfeature
	assert.is(x.fun, "function")
	arglis <- list(save.time = save.time, ...)
	portion.above.sum <- err.sum <- bias.sum <- score.sum <- as.numeric(rep(0, nfeature))
	P0.bias.sum <- ncp.bias.sum <- 0
	null.indicator.arg <- null.indicator
	reset.null.indicator <- length(null.indicator.arg) == 1 && is.character(null.indicator.arg)
	for(k in 1:nsim)
	{
		x <- x.fun()
		get.posteriorP0 <- function(fun, ...)
		{
			pos0 <- do.call(fun, c(arglis, list(x = x, ...)))
			assert.is(pos0, "posteriorP0")
			pos0
		}
		post0 <- get.posteriorP0(fun = posteriorP0.fun)
		true.post0 <- try(get.posteriorP0(fun = true.posteriorP0.fun, unknown.param = unknown.param))
		if(is.err(true.post0))
		{ message("bad true.post0"); browser()}
		if(!(is.prob(post0) && is.prob(true.post0)))
		{ message("bad get.posteriorP0 return value"); browser()}
		if(reset.null.indicator)
		{
			if(debug.null.indicator)
			{
				message("debugging null.indicator")
				if(null.indicator.arg == "x")
					print(data.frame(null.indicator = as(x@null.indicator, "numeric"), true.post0 = as(true.post0, "numeric")))
				browser()
			}
			null.indicator <- if(null.indicator.arg == "x")
			{
				assert.is(x, "r.twoDistrMixture") # required to have null.indicator slot
				x@null.indicator
			}
			else if(null.indicator.arg == "true.post0")
				ifelse(true.post0 >= 1 / 2, 1, 0)
			else
				stop("bad null.indicator.arg")
		}
		assert.is(null.indicator, "numeric")
		stopifnot(all(length(err.sum) == c(length(post0), length(true.post0), length(null.indicator))))
		discrepancy <- as(post0 - true.post0, "numeric")
		portion.above.sum <- portion.above.sum + ifelse(discrepancy >= 0, 1, ifelse(discrepancy == 0, 1/2, 0)) # >= on 110528
		err.sum <- err.sum + discrepancy ^ 2
		bias.sum <- bias.sum + discrepancy
		score.sum <- score.sum + (post0 - null.indicator) ^ 2
		if(!identical(unknown.param, unknownParam(true.post0))) # && !is.nothing(unknownParam(true.post0)))
		{
			message("bad unknown.param or unknownParam(true.post0)?"); print(unknown.param); print(unknownParam(true.post0))
			browser()
		}
		unknown.param.estimate <- unknownParam(post0)
		if(is.nothing(unknown.param.estimate)) # nonparametric density estimation
		{
			P0.bias.sum <- P0.bias.sum + NA
			ncp.bias.sum <- ncp.bias.sum + NA
		}
		else # parametric density estimation
		{
			bias.term <- function(name)
			{
				if(name %in% names(unknown.param))
				{
					stopifnot(name %in% names(unknown.param.estimate))
					term <- unknown.param.estimate[name][[1]] - unknown.param[name][[1]]
					stopifnot(length(term) == 1 && is.finite(term))
					term
				}
				else
					NA
			}
			if(!all(c("P0", "ncp") %in% names(unknown.param.estimate)))
			{
				message("bad names(unknown.param.estimate):")
				print(names(unknown.param.estimate))
				print(unknown.param.estimate)
				browser()
			}
			P0.discrepancy <- bias.term("P0")
			if(abs(P0.discrepancy) > max.P0.discrepancy || (unknown.param$P0 < 1 && abs(P0.discrepancy) < min.P0.discrepancy))
			{
				message("P0.discrepancy:", P0.discrepancy); print(unknown.param.estimate); print(unknown.param)
				message("consider running lfun <- likFun(object = post0@family, x = post0@x)")
				message("consider running rpost0=do.call(resampled.posteriorP0, c(list(x = x, nsample = 1), arglis))")
				browser()
			}
			P0.bias.sum <- P0.bias.sum + P0.discrepancy
			ncp.bias.sum <- ncp.bias.sum + bias.term("ncp")
		}
	}
	stopifnot(all(is.finite(err.sum)))
	names(portion.above.sum) <- names(err.sum) <- names(bias.sum) <- names(P0.bias.sum) <- names(ncp.bias.sum) <- names(score.sum) <- names(null.indicator) <- NULL
	sim.mean <- function(total){total / nsim}
	new.posteriorP0.error(sim.mean(err.sum), portion.above = sim.mean(portion.above.sum), bias = sim.mean(bias.sum), P0.bias = sim.mean(P0.bias.sum), ncp.bias = sim.mean(ncp.bias.sum), score = sim.mean(score.sum), null.indicator = null.indicator, unknown.param = unknown.param, ann = ann, estimator = posteriorP0.fun, x.fun = x.fun)
} # end posteriorP0.error

new.posteriorP0.errors <- function(object, name = names(object))
{
	stopifnot(length(name) %in% c(0, length(object)))
	names(object) <- name
	errs <- new("posteriorP0.errors", object)
	names(errs) <- name
	errs
}

new.posteriorP0.errors.lis <- function(object, nfeature)
{
	new("posteriorP0.errors.lis", object, nfeature = Numeric(nfeature))
}

new.posteriorP0.errorPair <- function(uncorrected, corrected)
{
	new("posteriorP0.errorPair", uncorrected = uncorrected, corrected = corrected)
}
posteriorP0.errorPair <- function(..., method = default("p", "posteriorP0.errorPair method"), nresample = NULL)
{
	arglis <- list(...)
	if("posteriorP0.fun" %in% names(arglis))
		stop("bad posteriorP0.errorPair args")
	message("\n\n\n\nCOMPUTING posteriorP0.errorPair")
	message("\n\nComputing uncorrected LFDR performance on ", date())
	uncorrected <- posteriorP0.error(posteriorP0.fun = posteriorP0, ...)
	message("\n\nComputing corrected LFDR performance on ", date())
	corrected <- posteriorP0.error(method = method, nresample = nresample, ...)
	new.posteriorP0.errorPair(uncorrected = uncorrected, corrected = corrected)
}


setMethod("P0", signature(object = "lik.ratio"), function(object, correction = 0, P0 = numeric(0), locationFUN, ...)
{
	stopifnot(is.prob(correction))
	if(is.nothing(P0) && length(list(...)) == 0 && missing(locationFUN)) # MLE
	{
		null.boo <- codelength(object) >= 0 # sign error corrected 100808
		if(any(is.na(null.boo)))
			stop("null.boo problem")
		Scalar((sum(null.boo) + correction) / (length(null.boo) + 2 * correction))
	}
	else if(is.prob(P0) && length(P0) == 1)
	{
		if(missing(locationFUN))
			locationFUN <- mean
		if(identical(locationFUN, "hard"))
			locationFUN <- function(x){mean(ifelse(x >= 0.5, 1, 0))}
		assert.is(locationFUN, "function")
		lfdr <- LFDR(object = object, P0 = P0, ...)
		stopifnot(is.prob(lfdr))
		location <- locationFUN(lfdr)
		if(!is.prob(location))
		{ message("unexpected location: ", location); browser()}
		Scalar(location)
	}
	else
		stop("bad P0 argument")
})
setMethod("P0", signature(object = "lik.ratios"), function(object, ...)
{
	sapply(object, P0, ...)
})
#-----------------
lik.ratio.odds <- function(object, base = exp(1))
{
	vec <- log(odds(object), base = base)
	new.lik.ratio(vec, type = "odds", base = base, unknown.param.set = list())
}
posterior.odds <- function(object, P0 = numeric(0), ...)
{
	assert.is(object, "lik.ratio")
	stopifnot(length(object@base) == 1 && object@base > 0)
	stopifnot(length(list(...)) == 0)
	if(is.nothing(P0))
	{
		P0 <- P0(object)
		P1 <- P1(object)
	}
	else if(is.prob(P0) && length(P0) == 1)
		P1 <- 1 - P0
	else
		stop("bad P0")
	stopifnot(P0 + P1 == 1)
	bf1 <- object
	stopifnot(length(bf1) == length(object))
	prior.odds1 <- log(P1 / P0, base = object@base)
	post.odds1 <- prior.odds1 + bf1
	names(post.odds1) <- names(object)
	assert.is(post.odds1, "lik.ratio")
	post.odds1@type = "posterior odds"
	post.odds1
}
LFDR <- function(object, ...)
{
	if(is(object, "lik.ratio"))
	{
		stopifnot(length(object@base) == 1 && object@base > 0)
		post.Odds0 <- object@base ^ (-posterior.odds(object = object, ...))
		vec <- as(1 / (1 + 1 / post.Odds0), "numeric")
		stopifnot(is.prob(vec))
		stopifnot(length(vec) == length(object))
		names(vec) <- names(object)
		stopifnot(is.prob(vec))
		Numeric(vec)
	}
	else if(is(object, "lik.ratios"))
		lapply(object, LFDR, ...)
	else if(is(object, "list"))
		LFDR(lik.ratios(object), ...)
	else
		stop("cannot compute the local false discovery rate")
}
lower.LFDR <- function(object, P0 = 1 / 2, ...) # moved from "combine111006.r" on 111025
{
	ev.name <- "ev.ban"
	if(is(object, "data.frame") && ev.name %in% names(object))
	{
		lr0 <- 1 / 10 ^ object[, ev.name]
		names(lr0) <- nam <- rownames(object)
	}
	else if(is(object, "lik.ratio"))
		return(LFDR(object = object, P0 = P0, ...))
	else
		stop("not yet implemented for that argument")
	assert.is(lr0, "numeric")
	na.boo <- is.na(lr0)
	if(any(na.boo))
	{
		message("removing ", sum(na.boo), " missing features")
		lr0 <- lr0[!na.boo]
	}
	lr0 <- Numeric(lr0)
	post.odds0 <- lr0 * P0 / (1 - P0)
	post.p0 <- post.odds0 / (1 + post.odds0)
	stopifnot(is.prob(post.p0))
	names(post.p0) <- nam[!na.boo]
	if(is.null(names(post.p0)))
	{ message("bad names(post.p0)"); browser()}
	post.p0
}
is.plausible <- function(object, probs, ...)
{
	if(!is(probs, "probs.to.combine")) # combine.r
		warning("probs is of unexpected class")
	stopifnot(is(probs, "list"))
	nam <- names(probs)
	lower.bound <- lower.LFDR(object, ...)[nam]
	stopifnot(sameNames(lower.bound, probs))
	vec <- sapply(probs, function(prob){all(prob[nam] >= lower.bound)})
	names(vec) <- probs@method.name
	vec
}
##
#
#
## near end of file:
#
##when.MDL.last.loaded <- date()
#
## EOF
