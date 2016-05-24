# NMWL.s created by David Bickel on 11 September 2010 by moving material from MDL.s.

approxWeight <- function(object, x)
{
	assert.is(object, "Weight")
	weight <- object
	w.incidental <- sum(weight@incidental.weight)
	incidental.weight <- scalar(weight@incidental.weight[1])
	stopifnot(all(weight@incidental.weight == incidental.weight))
	ncomp <- ncomparison(weight) # old number of comparisons is new number of incidental weights
	stopifnot(length(weight@incidental.weight) == ncomp - 1)
	incidental.weight <- scalar(incidental.weight * (ncomp - 1) / ncomp)
	weight@incidental.weight <- rep(incidental.weight, ncomp)
	stopifnot(ncomparison(weight) == ncomp + 1)
	stopifnot(ncomp == length(weight@incidental.weight))
	stopifnot(abs(sum(weight@incidental.weight) - w.incidental) < min(1e-3, 1 / ncomp))
	stopifnot(validObject(weight))
	if(!missing(x))
	{
		assert.is(x, "statistic")
		stopifnot(length(x) == length(weight@incidental.weight))
	}
	weight
}
Weight <- function(x, incidental, ncomparison, size = numeric(0), incidental.size = numeric(0), lower.bound = numeric(0), upper.bound = numeric(0)) # "size" is that of focus comparison
{
	if(!missing(ncomparison) && !is.nothing(size) && missing(x) && missing(incidental))
	{
		get.Weight <- function(focus.weight, incidental.weight)
		{
			new.Weight(focus.weight = focus.weight, incidental.weight = incidental.weight, lower.bound = lower.bound, upper.bound = upper.bound)
		}
		N.ok <- function(N){length(N) == 1 && N >= 1}
		stopifnot(N.ok(ncomparison) && ncomparison == floor(ncomparison))
		stopifnot(N.ok(size))
		w.incidental <- 1 / (size + 1)
		num.incidental <- ncomparison - 1
		incidental.weight <- if(is.nothing(incidental.size))
		{
			rep(w.incidental / num.incidental, num.incidental)
		}
		else if(num.incidental == length(incidental.size))
			w.incidental * incidental.size / sum(incidental.size)
		stopifnot(abs(sum(incidental.weight) - w.incidental) < min(1e-3, 1 / num.incidental))
		get.Weight(focus.weight = 1 - w.incidental, incidental.weight = incidental.weight)
	}
	else if(missing(ncomparison) && is.nothing(size) && !missing(x))
	{
		if(missing(incidental)) incidental <- FALSE
		if(is(x, "xprnSet"))
		{
			mat <- exprs(x)
			ncomparison <- nrow(mat)
			size <- ncol(mat)
		}
		else if(is(x, "xprnSetPair"))
		{
			xmat <- exprs(x@x)
			ymat <- exprs(x@y)
			ncomparison <- nrow(xmat)
			stopifnot(nrow(ymat) == ncomparison)
			size <- ncol(xmat) + ncol(ymat)
		}
		else
			stop("bad x class")
		if(incidental)
			ncomparison <- ncomparison + 1
		Weight(ncomparison = ncomparison, size = size, lower.bound = lower.bound, upper.bound = upper.bound)
	}
	else
		stop("cannot compute Weight from arguments given")
}


UNML <- function(object, mle.fun, truncation.tol = numeric(0), ...)
{
	if(missing(object))
	{
		if(missing(mle.fun))
			stop("not enough arguments") # missing(mle.fun) <- MLE
		assert.is(mle.fun, "function")
		Lik <- mle.fun(..., return.MLE = FALSE, log = FALSE, truncation.tol = truncation.tol)
		if(!is.finite(Lik))
		{ message("error in UNML"); browser()}
		Lik
	}
	else if(missing(mle.fun) && is(object, "extendedFamily"))
		UNML(mle.fun = object@mle.fun, truncation.tol = truncation.tol, ...)
	else if(missing(mle.fun) && is(object, "finiteFamily") && length(truncation.tol) == 0)
	{
		mle.fun <- function(...)
		{
			MLE(object = object, ...)
		}
		UNML(mle.fun = mle.fun, truncation.tol = numeric(0), ...) # mle.fun = object@mle.fun
	}
	else
		stop("bad UNML args")
} # end UNML
new.DensityFunUNML <- function(object, weight)
{
	new("DensityFunUNML", object, weight = weight)
}
DensityFunUNML <- function(object = stop("object missing"), truncation.tol = stop("truncation.tol missing"), incidental.x, weight = stop("weight missing")) # This low-level function as a rigid argument list by design; beware of making it flexible.
{
	if(missing(incidental.x))
	{ message("missing incidental.x"); browser()}
	if(!is(incidental.x, "numeric"))
	{ message("wrong incidental.x"); browser()}
	assert.is(incidental.x, "numeric")
	assert.is(weight, c("numeric", "Weight"))
	if(!is(weight, "Weight"))
		weight <- new.Weight(incidental.weight = weight)
	component.UNML <- function(component.x) # rigid argument list by design; beware of making it flexible
	{
		assert.is(component.x, "numeric")
		sapply(component.x, function(x){UNML(object = object, x = x, truncation.tol = truncation.tol, incidental.x = incidental.x, weight = weight)})
	}
	fun <- function(x) # rigid argument list by design; beware of making it flexible
	{
		assert.is(x, "numeric")
		sapply(x, function(dummy.x)
		{
			dummy.x.unml <- component.UNML(component.x = dummy.x)
			stopifnot(length(dummy.x.unml) == 1)
			dummy.x.unml
		})
	}
	new.DensityFunUNML(fun, weight = weight)
}

new.Complexity <- function(object, weight)
{
	new("Complexity", Scalar(object), weight = weight)
}
Complexity <- function(object, lower.x, upper.x, p = 1e-3, truncation.tol = numeric(0), allow.0.Complexity = FALSE, crude = FALSE, incidental.x = numeric(0), weight = numeric(0), ...)
{
	if(is(object, "DensityFunNML"))
		object@COMP
	else if(is(object, "DensityFunNMLs"))
	{
		vec <- sapply(object, Complexity)
		if(all(vec == vec[1]))
			Complexity(object[[1]])
		else
		{
			names(vec) <- names(object)
			vec
		}
	}
	else if(is(object, "Family"))
	{
		ok <- try(!object@discrete || is(object, "finiteFamily"))
		if(is.err(ok) || !ok)
		{	message("object must be continuous for numeric integration"); browser()}
		if(crude && is(object, "finiteFamily") && length(truncation.tol) == 0 && length(object@unknown.param.set) == 1 && length(object@unknown.param.set[[1]]) > 1 && missing(lower.x) && missing(upper.x) && is.nothing(weight))
		{
			stop("this crude approximation is extremely misleading unless there is sufficient separation, and it does not yet return Complexity object")
			subtot <- sapply(object@unknown.param.set[[1]], function(elem)
			{
				fam <- object
				fam@unknown.param.set <- list(elem)
				if(length(fam@unknown.param.set[[1]]) != 1)
				{ message("wrong length"); browser()}
				names(fam@unknown.param.set) <- names(object@unknown.param.set)
				stopifnot(validObject(fam))
				Complexity(object = fam, p = p, truncation.tol = numeric(0), allow.0.Complexity = allow.0.Complexity, ...)
			})
			assert.is(subtot, "numeric")
			stopifnot(length(subtot) == length(object@unknown.param.set[[1]]))
			sum(subtot)
		}
		else if(!crude && is(object, "extendedFamily") || (is(object, "finiteFamily") && length(truncation.tol) == 0 && length(object@unknown.param.set) == 1))
		{
			if("x" %in% names(list(...)))
				stop("x was specified for Complexity")
			if("param0" %in% names(list(...)))
				stop("param0 was specified for Complexity")
			is.ok <- function(x){length(x) == 0 || is.finite(x)}
			if(missing(lower.x) || missing(upper.x))
			{
				if(length(p) == 0)
				{
					fac <- 0.5
					if(missing(lower.x)) lower.x <- fac * object@lower.unknown.param
					if(missing(upper.x)) upper.x <- fac * object@upper.unknown.param
				}
				else if(is.ok(object@lower.unknown.param) && is.ok(object@upper.unknown.param) && length(p) == 1)
				{
					message("This could be a while. ", date())
					lims <- integration.limits(object = object, p = p) # slow
					if(missing(lower.x)) lower.x <- lims["lower"]
					if(missing(upper.x)) upper.x <- lims["upper"]
				}
				else
					stop("need lower.x and upper.x compute parametric complexity for the given family of distributions")
			}
			stopifnot(lower.x < upper.x)
			vector.UNML <- DensityFunUNML(object = object, truncation.tol = truncation.tol, incidental.x = incidental.x, weight = weight)
			integral <- try(integrate(f = vector.UNML, lower = lower.x, upper = upper.x, ...))
			if(is(integral, "try-error"))
			{ warning("bad integral"); message("bad integral"); browser()}
			if(integral$value == 0)
			{
				if(!allow.0.Complexity)
				{ message("parametric Complexity is 0"); browser()}
			}
			new.Complexity(integral$value, weight = as(vector.UNML, "Weight"))
		}
	} # end if(is(object, "Family")
	else
		stop("cannot determine the parametric complexity of the given family of distributions")
}
NML <- function(object, x, lower.x, upper.x, param0, COMP, truncation.tol = numeric(0), allow.0.Complexity = FALSE, incidental.x = numeric(0), weight = numeric(0), ...)
{
	assert.is(object, c("extendedFamily", "finiteFamily"))
	assert.is(x, "numeric")
	stopifnot(length(x) == 1) # for compatibility with the use of "integrate" in "Complexity" and "divergence"
	get.UNML <- DensityFunUNML(object = object, truncation.tol = truncation.tol, incidental.x = incidental.x, weight = weight)
	if(missing(COMP))
	{
		COMP <- Complexity(object = object, lower.x = lower.x, upper.x = upper.x, truncation.tol = truncation.tol, incidental.x = incidental.x, weight = weight, allow.0.Complexity = allow.0.Complexity, ...)
	}
	else if(!is(COMP, "Complexity") && is(COMP, "numeric") && missing(lower.x) && missing(upper.x))
		stop("numeric non-Complexity COMP is no longer accepted as an argument") # Scalar(COMP)
	else if(!is(COMP, "Complexity"))
		stop("NML cannot use the COMP provided as the parametric complexity of object")
	assert.is(COMP, "Complexity")
	numer <- get.UNML(x) / COMP
	if(missing(param0))
		numer
	else if(is.list(param0) || is.numeric(param0))
	{
		if(is.numeric(param0))
		{
			nam <- names(param0)
			param0 <- list(param0)
			names(param0) <- nam
		}
		assert.is(param0, "list")
		stopifnot(length(object@unknown.param.name) == length(object@unknown.param.name))
		if(is.null(names(param0)))
		{
			if(length(param0) == 1)
				names(param0) <- object@unknown.param.name
			else
				stop("bad param0")
		}
		stopifnot(setequal(names(param0), object@unknown.param.name))
		lf <- likFun(object = object, base = exp(1), x = x)
		lik0 <- do.call(lf, param0)
		Lik0 <- exp(lik0)
		Scalar(numer / Lik0)
	}
	else
		stop("bad NML args")
} # end NML

UNML.Binom <- function(x, size, approximate = FALSE)
{
	if(approximate)
	{
		fam <- extendedFamily.Binom(size = size)
		sapply(x, function(X){UNML(object = fam, x = X)})
	}
	else
	{
		assert.are(list(x, size), "numeric")
		stopifnot(length(size) == 1 && size >= 1)
		stopifnot(all(x >= 0 && x <= size))
		get.Lik <- function(prob)
		{
			Lik.Binom(x = x, size = size, prob = prob)
		}
		get.Lik(prob = x / size) #dbinom(x = x, size = size, prob = x / size, log = FALSE)
	}
}
Complexity.Binom <- function(size, min.prob = 0, max.prob = 1, call.entropy = FALSE)
{
	assert.is(size, "numeric")
	stopifnot(length(size) == 1 && size >= 1)
	stopifnot(scalar(min.prob) >= 0 && scalar(max.prob) <= 1 && min.prob <= max.prob)
	x <- 0:size
	prob.hat <- x / size
	boo <- prob.hat >= min.prob & prob.hat <= max.prob
	k <- x[boo]
	vec <- if(call.entropy && min.prob == 0 && max.prob == 1)
	{
		stop("error in this not corrected")
		prob <- prob.hat[boo]
		stopifnot(length(k) == length(prob))
		choose(n = size, k = k) * exp(- size * entropy(P = prob, base = exp(1))) # Grünwald, 2007, prob. 227
	}
	else
		UNML.Binom(x = k, size = size)
	stopifnot(length(vec) == sum(boo))
	COMP <- sum(vec)
	if(is.na(COMP))
	{ message("bad COMP"); browser()}
	COMP
}
NML.Binom <- function(x, size, min.prob = 0, max.prob = 1, prob0, plot = missing(x) && length(size) == 1, base = 2, ...)
{
	if(missing(x) && length(size) == 1)
		x <- 0:size
	assert.are(list(x, size), "numeric")
	stopifnot(length(size) == 1 && size >= 1)
	stopifnot(scalar(min.prob) >= 0 && scalar(max.prob) <= 1 && min.prob <= max.prob)
	prob.hat <- x / size
	boo <- prob.hat >= min.prob & prob.hat <= max.prob
	numerator <- ifelse(boo, UNML.Binom(x = x, size = size), 0) / Complexity.Binom(size = size, min.prob = min.prob, max.prob = max.prob, ...)
	Nml <- if(missing(prob0))
		numerator
	else if(is.prob(prob0) && (length(numerator) == 1 || length(prob0) == 1 || length(numerator) == length(prob0)))
	{
		get.Lik <- function(prob)
		{
			Lik.Binom(x = x, size = size, prob = prob)
		}
		numerator / get.Lik(prob = prob0)
	}
	else
		stop("bad prob0")
	if(plot)
	{
		ylab <- if(base == 2)
			"bits"
		else if(base == 10)
			"bans"
		else
			"log"
		plot(x = x, y = logb(Nml, base = base), log = "", ylab = ylab)
	}
	Nml
}

new.DensityFunNML <- function(object, family, COMP) #, weight = blank.Weight())
{
	new("DensityFunNML", object, family = family, COMP = COMP) # , weight = weight)
}
DensityFunNML <- function(object, COMP = NULL, truncation.tol = numeric(0), incidental.x = numeric(0), weight = numeric(0), allow.0.Complexity = FALSE, verbose, allow.unused.args = FALSE, ...)
{
	assert.is(object, c("extendedFamily", "finiteFamily"))
	if("x" %in% names(list(...)))
		stop("x was specified for DensityFunNML")
	if("param0" %in% names(list(...)))
		stop("param0 was specified for DensityFunNML")
	if(is(incidental.x, "xprnSetObject") && is.nothing(weight)) #  && !missing(i) && length(i) == 1)
	{
		weight <- Weight(x = incidental.x, incidental = TRUE) # as(incidental.x, "Weight")
		es <- incidental.x # [-i, ]
		#XXX: from statistic to nstatistic
		incidental.x <- nstatistic(x = es, object = object)
		stopifnot(length(incidental.x) == length(weight@incidental.weight))
	}
#	else
#		stopifnot(missing(i))
	assert.is(incidental.x, "numeric")
	assert.is(weight, c("numeric", "Weight"))
	COMP <- if(is.nothing(COMP))
	{
		if(missing(verbose)) verbose <- default(TRUE, "verbose")
		Complexity(object = object, truncation.tol = truncation.tol, incidental.x = incidental.x, weight = weight, allow.0.Complexity = allow.0.Complexity, ...)
	}
	else if(length(list(...)) == 0 || allow.unused.args)
	{
		if(missing(verbose)) verbose <- FALSE
		if(is(COMP, "Complexity"))
			COMP
		else if(is(COMP, "numeric"))
			stop("numeric non-Complexity COMP is not accepted as of 100916") # Scalar(COMP)
		else
			stop("non-numeric COMP")
	}
	else
		stop("DensityFunNML cannot use the COMP provided as the parametric complexity of object")
	if(verbose)
		message("COMP: ", COMP)
	assert.is(COMP, "Complexity", "outside DensityFun")
	DensityFun <- function(x, ...)
	{
		assert.is(COMP, "Complexity", "inside DensityFun")
		NML(object = object, x = x, COMP = COMP, truncation.tol = truncation.tol, incidental.x = incidental.x, weight = weight, allow.0.Complexity = allow.0.Complexity, ...)
	}
	new.DensityFunNML(object = DensityFun, family = object, COMP = COMP)
} # end DensityFunNML

ncp1bdd.t1 <- 1.5666625355574042 # from "noncentral distributions.nb"
distinct.unknown.param <- function(object, lower, upper, param0, ...)
{
	stopifnot(is(object, "Family") && length(object@unknown.param.name) == 1 && length(param0) == 1)
	param.value0 <- if(is(param0, "list"))
	{
		if(length(param0) == 1 && sameSet(names(param0), object@unknown.param.name))
			param0[[1]]
		else
			stop("bad parameter list")
	}
	else
		param0
	assert.is(param.value0, "numeric")
	stopifnot(length(param.value0) == 1)
	get.unknown.param <- function(param.value)
	{
		assert.is(param.value, "numeric")
		stopifnot(length(param.value) == 1)
		unknown.param <- list(param.value)
		names(unknown.param) <- object@unknown.param.name
		unknown.param
	}
	get.distr <- function(param.value)
	{
		unknown.param <- get.unknown.param(param.value)
		Distr(object = object, unknown.param = unknown.param)
	}
	f <- function(param.value1)
	{
		distr1 <- get.distr(param.value = param.value1)
		distr0 <- get.distr(param.value = param.value0)
	# Quantile[NoncentralStudentTDistribution[\[Nu], \[Delta]], .25] - Quantile[NoncentralStudentTDistribution[\[Nu], 0], .75] (*noncentralStudentTBDD of "noncentral distributions.nb"*)
		q(distr1)(0.25) - q(distr0)(0.75)
	}
	get.unknown.param(uniroot(f = f, lower = lower, upper = upper, ...)$root)
}
distinct.unknown.param.Td <- function(df, ncp0 = 0)
{
	distinct.unknown.param(object = Family.Td(df), lower = 1e-2 + ncp0, upper = 2 + ncp0, param0 = list(ncp = ncp0))
}
DensityFunNML.Td <- function(size, df, local.alternative, max.ncp, ncp0 = 0, ncp.set.i = 1:2, ...)
{
	if(missing(df))
		df <- size - 1
	else
		stopifnot(missing(size))
	if(missing(max.ncp))
	{
		max.ncp <- if(local.alternative)
			distinct.unknown.param.Td(df = df, ncp0 = ncp0)$ncp
		else if(ncp0 == 0)
			ncp1bdd.t1 * sqrt(size)
		else
			stop("incompatible arguments")
	}
	else if(!missing(local.alternative))
		stop("cannot determine noncentrality parameter from the arguments given")
	assert.is(max.ncp, "numeric")
	max.ncp <- Scalar(max.ncp)
	fun <- DensityFunNML(finiteFamily.Td(df = df, ncp.set = c(-max.ncp, max.ncp)[ncp.set.i]), ...)
	fun0 <- as(fun, "function")
	fun@.Data <- function(x, x.converges = FALSE, ...)
	{#XXX: from statistic to nstatistic
		x <- nstatistic(object = fun@family, x = x, x.converges = x.converges, size = size)
		fun0(x = x, ...)
	}
	fun
}
optimize.Td <- function(...)
{
	df <- 1
	get.dis <- function(max.ncp)
	{
		P1 <- DensityFunNML(extendedFamily.Td(df = df, lower.ncp = -max.ncp, upper.ncp = max.ncp))
		divergence(P1 = P1, P2 = Td(df = df, ncp = 0))
	}
	optimize(f = get.dis, ...)
}

new.DensityFunNMLs <- function(object, statistic.fun, reduced.x, weight, family)
{
	new("DensityFunNMLs", object, statistic.fun = statistic.fun, reduced.x = new.statistic(reduced.x), weight = weight, family = family)
}
DensityFunNMLs <- function(x, incidental.x = numeric(0), object, div.mult.factor, weight, size = numeric(0), statistic.fun, lower.x, upper.x, verbose.i, verbose = FALSE, save.time = default(FALSE, "save.time"), COMP, ...)
{
	if(missing(object))
	{
		object <- if(is(x, "xprnSetObject"))
		{
			if(missing(div.mult.factor))
				Family(x)
			else
				Family(x, div.mult.factor = div.mult.factor) # added 101129
		}
		else
			stop("cannot generate object")
	}
	else
		stopifnot(missing(div.mult.factor))
	assert.is(object, c("Family", "Families"))
	if(verbose && !missing(weight)) {message("0 weight:"); print(try(weight))}
	if(missing(weight))
	{
		get.Weight <- function(focus.size, incidental.size = numeric(0))
		{
			focus.size <- Scalar(focus.size)
			w <- try(if(is(x, "numeric") && length(focus.size) == 1)
				Weight(ncomparison = length(x), size = focus.size, incidental.size = incidental.size)
			else if(is(x, "xprnSetObject") && is.nothing(focus.size) && is.nothing(incidental.size))
				Weight(x = x)
			else
				stop("cannot compute weight"))
			if(is.err(w) || is.err(try(validObject(w))))
			{ message("try 1 failed"); browser()}
			if(verbose) {message("w:"); print(w)}
			w
		}
		assert.is(size, "numeric")
		weight <- if(is.nothing(size) || length(size) == 1)
			get.Weight(focus.size = size)
		else if(length(size) > 1)
			new.Weights(lapply(1:length(size), function(i)
			{
				get.Weight(focus.size = size[i], incidental.size = size[-i])
			}))
		else
			stop("the impossible ocurred")
		if(verbose) {message("1 weight:"); print(weight)}
	}
	else
		stopifnot(is.nothing(size))
	tri <- try(assert.is(weight, c("Weight", "Weights")))
	if(is.err(tri))
	{ message("try 2 failed"); browser()}
	ncomp <- ncomparison(weight) # 1 + length(weight@incidental.weight)
	stopifnot(ncomp >= 2)
	if(missing(statistic.fun))
	{
		statistic.fun <- if(is(x, "numeric"))
			function(x){new.statistic(x)} # no data reduction if reduced data are given
		else if(is(x, "xprnSetObject"))
		{
			stopifnot(ncomp == ncomparison(x))
			function(x){nstatistic(x = x, object = object)}
			#XXX: from statistic to nstatistic
		}
		else
			stop("cannot compute reduced data")
	}
	assert.is(statistic.fun, "function")
	reduced.x <- statistic.fun(x)
	rm(x)
	tri <- try(assert.is(reduced.x, "statistic"))
	if(is.err(tri))
	{ message("reduced.x failed"); browser()}
	stopifnot(ncomp == ncomparison(reduced.x)) # stopifnot(length(reduced.x) == ncomp)
	if(missing(lower.x) || missing(upper.x))
	{
		message("reduced.x:")
		print(stats(reduced.x))
		stop("lower.x or upper.x not specified")
	}
	get.incidental.x <- function(alternative.incidental.x)
	{
		ix <- if(is.nothing(incidental.x))
			alternative.incidental.x
		else
		{
			if(length(incidental.x) == 1)
				rep(incidental.x, length(alternative.incidental.x))
			else
				incidental.x
		}
		stopifnot(length(ix) == length(alternative.incidental.x))
		ix
	}
	if(missing(save.time))
		save.time <- default(is.nothing(incidental.x), "save.time") # default(FALSE, "save.time") changed 101119 # !is.nothing(incidental.x)
	if(missing(COMP))
	{
		COMP <- if(save.time)
		{
			if(is(weight, "Weights"))
				stop("cannot save time with comparison-dependent weights")
			assert.is(weight, "Weight")
			if(is(object, "Families"))
				stop("cannot save time with comparison-dependent families")
			assert.is(object, "Family")
			COMP <- Complexity(object = object, incidental.x = get.incidental.x(reduced.x), weight = approxWeight(weight, x = reduced.x), lower.x = lower.x, upper.x = upper.x, ...)
			message("COMP: ", COMP)
			COMP
		}
		else
			NULL
	}
	else
		stopifnot(is.nothing(incidental.x))
	if(missing(verbose.i))
	{
		verbose.i <- if(is.nothing(COMP))
			c(1:20, 1e2, 1e3, 1e4, 1e5)
		else
			numeric(0)
	}
	lis <- lapply(1:ncomp, function(i)
	{
		if(i %in% verbose.i)
			message("\nbeginning iteration ", i, " on ", date())
		incidental.x <- reduced.x[-i]
		object1 <- if(is(object, "Family"))
			object
		else if(is(object, "Families") && length(object) == ncomp)
			object[[i]]
		else
			stop("no object1 for DensityFunNML")
		weight1 <- if(is(weight, "Weight"))
			weight
		else if(is(weight, "Weights") && length(weight) == ncomp)
			weight[[i]]
		else
			stop("no weight1 for DensityFunNML")
		DensityFunNML(object = object1, incidental.x = get.incidental.x(incidental.x), weight = weight1, lower.x = lower.x, upper.x = upper.x, COMP = COMP, allow.unused.args = TRUE, ...)
	})
	names(lis) <- names(reduced.x)
	object1 <- if(is(object, "Family"))
		object
	else if(is(object, "Families") && length(object) == ncomp)
		blank.Family()
	else
		stop("no object1 for DensityFunNMLs")
	weight1 <- if(is(weight, "Weight"))
		weight
	else if(is(weight, "Weights") && length(weight) == ncomp)
		blank.Weight()
	else
		stop("no weight1 for DensityFunNMLs")
	new.DensityFunNMLs(lis, statistic.fun = statistic.fun, reduced.x = reduced.x, weight = weight1, family = object1)
}
es.DensityFunNMLs <- function(es, upper.x, ...) # from nml100916.r
{
	assert.is(es, "xprnSetObject")
	DensityFunNMLs(x = es, lower.x = 0, upper.x = upper.x, ...)
}

# near end of file:

when.NMWL.last.loaded <- date()

# EOF
