# interval.s created by David Bickel on 24 December 2010.

new.posteriorInterval <- function(lower, upper)
{
	stopifnot(length(lower) == length(upper))
	if(!(is.nothing(lower) && is.nothing(upper)))
	{
		if(any(is.na(c(lower, upper))))
			stop("missing values in confidence intervals")
		if(is.null(names(lower)))
			names(lower) <- names(upper)
		else if(is.null(names(upper)))
			names(upper) <- names(lower)
	}
	new("posteriorInterval", lower = lower, upper = upper)
}
#t.posteriorInterval <- function(p1 = 1 / 40, p2 = 1 - p1)#XXX|:no visible binding for global variable ÔlowerÕ and 'upper'
#{
#	stop("make CDF and then call p1 = 1 / 40, p2 = 1 - p1")
#	ci <- new.posteriorInterval(lower = lower, upper = upper)
#	stopifnot(is(ci, "posteriorInterval"))
#	ci
#}
posteriorInterval.apply <- function(...)
{
	mat <- sapply(...)
	assert.is(mat, "matrix")
	as(mat, "posteriorInterval")
}
posteriorInterval <- function(object, p1 = 1 / 40, p2 = 1 - p1)
{
	if(is(object, "inverseCDF"))
	{
		# posteriorInterval(as(object, "inverseCDFs"), p1 = p1, p2 = p2)
		ok <- length(p1) == 1 && length(p2) == 1 && p1 <= p2
		if(!ok)
		{ message("problem with p1 and p2"); browser()}
		int <- try(new.posteriorInterval(lower = object(p1), upper = object(p2)))
		if(is.err(int))
		{
			warning("bad inverseCDF object; returning NAs") # { message("bad inverseCDF object?"); browser()}
			new.posteriorInterval(lower = numeric(0), upper = numeric(0))
		}
		else
			int
	}
	else if(is(object, "inverseCDFs"))
	{
		boo <- !is.nothing(object)
		stopifnot(length(boo) == length(object))
		stopifnot(any(boo))
		if(!all(boo))
		{
			message("\n\nposteriors present only for ", sum(boo), " of ", length(boo), " attempted posteriors:\n")
			object <- object[boo]
			assert.is(object, "inverseCDFs")
		}
		stopifnot(length(object) == sum(boo))
		lis <- lapply(object, posteriorInterval, p1 = p1, p2 = p2)
		names(lis) <- names(object)
		if(!all(sapply(lis, is, class2 = "posteriorInterval")))
		{ message("wrong class"); browser()}
		boo <- !sapply(lis, is.nothing)
		stopifnot(length(boo) == length(lis))
		stopifnot(any(boo))
		if(!all(boo))
		{
			message("\n\nposterior intervals computed for ", sum(boo), " of ", length(boo), " posteriors:\n")
			lis <- lis[boo]
			assert.is(lis, "list")
		}
		as(lis, "posteriorInterval")
	}
	else
		stop("cannot compute posteriorInterval")
}

progress.index <- function(object, nreport = floor(min(30, len)))
{
	len <- if(length(object) == 1 && is(object, "numeric"))
		object
	else
		length(object)
	progress.i <- floor(sort(unique(c(1:5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 1e4, 2e4, 5e4, seq(min(10, len), len, length.out = nreport)))))
	progress.i[progress.i %in% 1:len]
}

Mass0 <- function(object)
{
	lfdr(object, max.lfdr = 1)
}

new.posteriorIntervalPlus <- function(object, alternative, param, p1, p2)
{
	add.nam <- function(y)
	{
		if(is.null(names(y)))
			names(y) <- names(object)
		y
	}
	new("posteriorIntervalPlus", object, alternative = add.nam(alternative), param = add.nam(param), p1 = Scalar(p1), p2 = Scalar(p2))
}
posteriorIntervalPlus <- function(size, nfeature, P0, nsample, p1 = 1 / 40, p2 = 1 - p1, mass0.fun = assumedNull, nfeatureIID = nfeature, ...)
{
#	if(length(p1) == 1) p1 <- rep(p1, nfeature)
#	if(length(p2) == 1) p2 <- rep(p2, nfeature)
	p1 <- Scalar(p1)
	p2 <- Scalar(p2)
	alternative <- as.logical(rep(NA, nsample))
	param <- as.numeric(rep(NA, nsample))
#	int.mat <- matrix(rep(NA, 2 * nfeature), nrow = 2)
	int.lis <- list()
	progress.i <- progress.index(nsample)
	for(i in 1:nsample)
	{
		if(i %in% progress.i)	message("\n\nstarting simulation number ", i, " of ", nsample, " on ", date())
		es <- Bernoulli.normal.rxprnSet(size = size, nfeature = nfeature, P0 = P0, nfeatureIID = nfeatureIID)
		if(i %in% progress.i)	message("  fixed-parameter CI ", i, " of ", nsample, " on ", date())
		j <- 1 # feature index
		icdfs.alt <- inverseCDFs(es[j, ], ...)
		if(i %in% progress.i)	message("  computing hypothesis probability vector ", i, " of ", nsample, " on ", date())
		mass0 <- if(is(mass0.fun, "function"))
			mass0.fun(es)
		else if(is(mass0.fun, "numeric"))
		{
			if(length(mass0.fun) == 1)
				mass0.fun <- rep(mass0.fun, nfeature)
			stopifnot(is.prob(mass0.fun) && length(mass0.fun) == nfeature)
			mass0.fun
		}
		if(i %in% progress.i)	message("  adding hypothesis probability vector ", i, " of ", nsample, " on ", date())
		icdf <- addPointMass(icdfs.alt[[j]], mass0 = if(is.prob(mass0)) mass0[j] else Mass0(mass0)[j])
		if(i %in% progress.i)	message("  random-parameter CI ", i, " of ", nsample, " on ", date())
		int <- posteriorInterval(icdf, p1 = p1, p2 = p2)
		int.lis <- c(int.lis, list(int))
		alternative[i] <- es@alternative[j]
		param[i] <- es@param[j]
	}
	len.ok <- length(alternative) == nsample && length(param) == nsample
	ok <- all(is.finite(alternative) & is.finite(param)) && length(int.lis) == nsample && len.ok
	if(!ok)
	{ message("bad posteriorIntervalPlus sims"); browser()}
	names(int.lis) <- make.names(paste("realization", 1:nsample, sep = "."))
	int <- as(int.lis, "posteriorInterval")
	new.posteriorIntervalPlus(int, alternative = alternative, param = param, p1 = p1, p2 = p2)
}


Width <- function(object, base = default(2, "base"))#change case by marta Nov2013
{
	assert.is(object, "posteriorInterval")
	vec <- object@upper - object@lower
	stopifnot(length(object) == length(vec))
	names(vec) <- names(object)
	assert.is(vec, "numeric")
	stopifnot(all(vec >= -1e-3))
	if(is.nothing(base))
		vec
	else
		new.log(vec, base = base)
}

covers <- function(object)
{
	assert.is(object, "posteriorIntervalPlus")
	vec <- object@lower <= object@param & object@param <= object@upper
	stopifnot(length(object) == length(vec))
	names(vec) <- names(object)
	assert.is(vec, "logical")
	vec
}
new.confidencePerformance <- function(coverage, coverage0, coverage1)
{
	add.name <- function(cover)
	{
		if(is.null(names(cover)))
			names(cover) <- names(coverage)
		Numeric(cover)
	}
	new("confidencePerformance", coverage = Numeric(coverage), coverage0 = add.name(coverage0), coverage1 = add.name(coverage1))
}
confidencePerformance <- function(object)
{
	if(is(object, "posteriorIntervalPlus"))
	{
		get.coverage <- function(pip)
		{
			sum(covers(pip)) / length(pip)
		}
		coverage <- get.coverage(object)
		coverage0 <- get.coverage(object[!object@alternative])
		coverage1 <- get.coverage(object[object@alternative])
	}
	else if(is(object, "posteriorIntervalPluses"))
	{
		stop("posteriorIntervalPluses not yet implemented")
		names(coverage) <- names(object)
	}
	else
		stop("cannot compute performance of confidence intervals")
	new.confidencePerformance(coverage = coverage, coverage0 = coverage0, coverage1 = coverage1)
}

new.inverseCDF <- function(object, min_param, max_param, param.name, type)
{
	new("inverseCDF", object, min_param = scalar(min_param), max_param = scalar(max_param), param.name = param.name, type = type)
}
#--------------

probability.inverseCDF<-function(...){}#XXX|:define this!!
#--------------
inverseCDF <- function(type, ...)
{
	fun <- if(type == "confidence")
		confidence.inverseCDF
	else if(type == "probability")
		probability.inverseCDF#XXX|:no visible binding for global variable Ôprobability.inverseCDFÕ
	else
		stop("bad type of inverseCDF")
	fun(...)
}
blank.inverseCDF <- function(object, param.name = "no parameter")
{
	if(missing(object))
		object <- "no function specified for inverseCDF"
	assert.is(object, "character")
	fun <- function(param)
	{
		stop(object)
	}
	new.inverseCDF(object = fun, min_param = -Inf, max_param = Inf, param.name = param.name, type = "no type")
}
confidence.inverseCDF <- function(object, conf.limit.fun, param.name, min_param, max_param, ...)
{
	if(missing(conf.limit.fun) && missing(param.name) && missing(min_param) && missing(max_param))
	{
		message("Due to missing arguments, calling t_test_inverseCDF.")
		t_test_inverseCDF(object = object, ...)
	}
	else
	{
		stopifnot(length(min_param) == 1 && length(max_param) == 1)
		if(missing(object))
		{ message("confidence.inverseCDF missing object"); browser()}
		inverse.cdf <- function(p) # integrate needs vector p
		{
			stopifnot(all(p >= 0 & p <= 1))
			sapply(p, conf.limit.fun, object = object)
		}
		new.inverseCDF(object = inverse.cdf, min_param = min_param, max_param = max_param, param.name = param.name, type = "confidence")
	}
} # end confidence.inverseCDF
t_test_inverseCDF <- function(object, ...)
{
	nsample <- if(is(object, "list") && length(object) == 2)
		2
	else if(is(object, "numeric"))
		1
	else
		stop(paste("Cannot assign nsample from given object of class", class(object)))
	sufficient.data <- if(is(object, "list") && length(object) == 2)
		Size(object[[1]]) >= 1 && Size(object[[2]]) >= 1 && Size(object[[1]]) + Size(object[[2]]) >= 3
	else if(is(object, "numeric"))
		Size(object) >= 2
	else
		stop(paste("Cannot assign sufficient.data from given object of class", class(object)))
	if(!sufficient.data)
		blank.inverseCDF()
	else
	{
		param.name <- if(nsample == 2)
			"difference of means"
		else if(nsample == 1)
			"mean"
		else
			stop(paste("cannot assign param.name from given object of class", class(object)))
		conf.limit.fun <- function(object, p)
		{
			get.conf.limit <- function(...)
			{
				t.test(..., conf.level = p, mu = 0, alternative = "less")$conf.int[2]
			}
			if(nsample == 2)
				get.conf.limit(x = object[[1]], y = object[[2]])
			else if(nsample == 1)
				get.conf.limit(x = object)
			else
				stop(paste("cannot compute", param.name, "p-value from given object of class", class(object)))
		}
		confidence.inverseCDF(object = object, conf.limit.fun = conf.limit.fun, param.name = param.name, min_param = -Inf, max_param = Inf, ...)
	}
}
z.test.inverseCDF <- function(object, prior.mean = numeric(0), complete = TRUE, mass0 = numeric(0), param0 = 0, ...)
{
	nsample <- if(is(object, "list") && length(object) == 2)
		2
	else if(is(object, "numeric"))
		1
	else
		stop(paste("Cannot assign nsample from given object of class", class(object)))
	sufficient.data <- if(is(object, "list") && length(object) == 2)
		stop("not implemented for two samples")
	else if(is(object, "numeric"))
		Size(object) >= 1
	else
		stop(paste("insufficient information from object of class", class(object)))
	if(!sufficient.data)
		blank.inverseCDF()
	else if(Size(object) == 1)
	{
		param.name <- "mean"
		marginal.confidence.interval <- complete && is.nothing(prior.mean)
		must.add.mass0 <- !is.nothing(mass0) && complete
		Mass0 <- if(must.add.mass0)
			mass0
		else
			numeric(0)
		conf.limit.fun <- function(object, p)
		{
			if(marginal.confidence.interval) # marginal confidence posterior
			{
				if(!is.nothing(mass0))
					stopifnot(must.add.mass0)
				qnorm(p = p, mean = object)
			}
			else
			{
				if(TRUE)
				{
					if(is.nothing(mass0) || must.add.mass0)
						mass0 <- 0
					else
						stopifnot(!complete) # !complete means not including probability mass at null
#					if(must.add.mass0)
#						stop("already incorporating mass")
					assert.is(mass0, "numeric")
					stopifnot(length(mass0) == 1)
					p1 <- ifelse(p < 1 / 2, p, 1 - p)
					alpha <- 2 * p1
					if(mass0 > 1 - alpha)
						param0
					else
					{
						rescale <- function(p1)
						{
							p1new <- p1 / (1 - mass0)
							if(!is.prob(p1new))
							{ print(c(p1 = p1, mass0 = mass0, p1new = p1new)); browser()}
							p1new
						}
						cond.limit <- function(p) # conditional on the alternative hypothesis
						{
							if(is.nothing(prior.mean)) # discontinuous confidence set
								qnorm(p = p, mean = object)
							else if(is.numeric(prior.mean)) # Efron (2010, p. 232); h110203
								qnorm(p = p, mean = (object + prior.mean) / 2, sd = 1 / sqrt(2))
							else
								stop("bad prior.mean")
						}
						cond.limit(p = ifelse(p < 1 / 2, rescale(p1 = p), 1 - rescale(p1 = 1 - p)))
					}
				}
			}
#			else
#				stop("bad z.test.inverseCDF arguments")
		}
		cicdf <- confidence.inverseCDF(object = object, conf.limit.fun = conf.limit.fun, param.name = param.name, min_param = -Inf, max_param = Inf, ...)
		if(must.add.mass0)
			addPointMass(object = cicdf, mass0 = mass0)
		else
			cicdf
	}
	else
		stop("object is not a single value")
}

plot_CI_nverseCDF <- function(object, mass0, p1 = 1 / 40, p2 = 1 - p1, prior.mean = numeric(0), complete = TRUE, param0 = 0, col = "gray", ...)
{
	assert.is(mass0, "numeric")
	if(is.nothing(mass0))
		stop("bad mass0")
	get.fun <- function(p)
	{
		function(mass0)
		{
			icdf <- if(is(object, "inverseCDF"))
				addPointMass(object = object, mass0 = mass0)
			else if(is(object, "numeric") && length(object) == 1)
			{
#				mass0 <- numeric(0)
				z.test.inverseCDF(object = object, prior.mean = prior.mean, complete = complete, mass0 = mass0, param0 = param0)
			}
			else
				stop("cannot get fun")
			icdf(p)
		}
	}
	rho <- p2 - p1
	ylab <- paste(percent(rho, return.numeric = FALSE), "posterior set")
	x <- if(is(object, "inverseCDF") || (is.nothing(prior.mean) && complete))
		mass0
	else
		mass0[mass0 <= rho]
	plot.polygon(x = x, y = list(get.fun(p = p1), get.fun(p = p2)), xlab = "local false discovery rate", ylab = ylab, col = col, xlim = range(mass0), in.ylim = param0, ...)
	if(!complete)
	{
		plot.polygon(x = mass0, y = rep(param0, length(mass0)), col = col, add = TRUE) #	abline(h = param0, col = col, cex = 2)
	}
	conditional.icdf <- z.test.inverseCDF(object = object, complete = TRUE, param0 = param0)
	abline(h = c(conditional.icdf(p1), conditional.icdf(p2)), lty = "dashed")
#	plot.polygon(x = c(0, 0), y = c(conditional.icdf(p1), conditional.icdf(p2)), col = "black", add = TRUE)
}

new.inverseCDFs <- function(object, size = numeric(0), mass0 = numeric(0))
{
	if(is(object, "inverseCDF"))
		object <- list(object)
	assert.is(object, "list")
	lis <- try(new("inverseCDFs", object, size = size, mass0 = Numeric(mass0)))
	if(is.err(lis) || !validObject(lis))
	{ message("bad new.inverseCDFs"); browser()}
	lis
}
inverseCDFs <- function(object, FUN, ...)
{
	if(is(object, "list") && missing(FUN) && is.nothing(list(...)))
		new.inverseCDFs(object)
	else if(is(object, "xprnSetObject"))
	{
		if(missing(FUN))
			FUN <- t_test_inverseCDF
		fun <- function(x){FUN(object = x, ...)}
		size <- Size(object)
		lis <- lapply(1:Length(object), function(i)
		{
			get.vec <- function(es)
			{
				assert.is(es, "xprnSet")
				vec <- as(exprs(es)[i, ], "numeric")
				if(is(es, "XprnSet")) logb(vec) else vec
			}
			x <- if(is(object, "xprnSet"))
				get.vec(object)
			else if(is(object, "xprnSetPair"))
				list(x = get.vec(object@x), y = get.vec(object@y))
			else
				stop("unrecognized expression object 1")
			ok <- is(x, "list") || is(x, "numeric")
			fun(x = x)
		})
		names(lis) <- names(object)
		stopifnot(sameNames(lis, size))
		new.inverseCDFs(object = lis, size = size)
	}
	else
		stop("inverseCDFs does not recognize your arguments")
}

setMethod("is.nothing", signature(object = "inverseCDF"), function(object)
{
	object@param.name == "no parameter" || object@type == "no type" || is.nothing(object@.Data)#marta added: || is.nothing(object@ .Data) Oct 2013
})
setMethod("is.nothing", signature(object = "inverseCDFs"), function(object)
{if(is.nothing(object@.Data)){return(TRUE)}#marta added: 
	boo <- sapply(object, is.nothing)
	names(boo) <- names(object)
	boo 
})
setMethod("is.nothing", signature(object = "posteriorInterval"), function(object)
{
	is.nothing(object@lower) && is.nothing(object@upper)
})

addPointMass <- function(object, mass0, param0 = 0, call.browser = FALSE, ...)
{
	if(is(object, "inverseCDF"))
	{
		assert.is(mass0, "numeric")
		icdf <- function(p) # integrate needs vector param
		{
			if(call.browser)
			{ message("calling the browser 1, as requested by the parameter call.browser"); browser()}
			mass0 <- Scalar(mass0)
			mass1 <- Scalar(1 - mass0)
			param0 <- scalar(param0)
			p <- Scalar(p)
			alpha <- Scalar(1 - p)
			param.low <- if(p <= mass1) object(p / mass1) else numeric(0)
			param.high <- if(alpha <= mass1) object(1 - alpha / mass1) else numeric(0)
#			stopifnot(param.low <= param.high)
			param <- if(length(param.low) == 1 && param.low < param0) # h101218f
				param.low
			else if(length(param.high) == 1 && param.high > param0)
				param.high
			else if(TRUE) # length(param.low) == 1 || length(param.high) == 1 # length(param.low) == 0 || length(param.high) == 0
				param0
			else
			{ message("that possibility was never considered"); browser()}
			if(is(param, "try-error"))
			{ message("What was wrong?"); browser()}
			if(call.browser)
			{ message("calling the browser 2, as requested by the parameter call.browser"); browser()}
			param
		}
		object2 <- object
		object2@.Data <- icdf
		if(!validObject(object2))
		{ message("bad ", class(object2)); browser()}
		object2
	}
	else if(is(object, "inverseCDFs"))
	{
		if(is(mass0, "empiricalNull"))
		{
			message("\ncomputing local false discovery rate estimates on ", date())
			mass0 <- Mass0(mass0)
		}
		assert.is(mass0, "numeric")
		if(length(mass0) == 1)
		{
			mass0 <- rep(mass0, length(object))
			names(mass0) <- names(object)
		}
		nam <- intersect(names(mass0), names(object))
		object2 <- object[nam]
		if(!is(object2, "inverseCDFs") || is.nothing(object2))
		{ message("object2 corrupt"); browser()}
#yo comente 18/10/2012		message("adding point masses on ", date())
		nreport <- floor(min(30, length(nam)))
#		progress.i <- floor(sort(unique(c(1:5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 1e4, 2e4, 5e4, seq(min(10, length(nam)), length(nam), length.out = nreports)))))
		progress.i <- progress.index(nam) #progress.i[progress.i %in% 1:length(nam)]
		progress.nam <- nam[progress.i]
		lis <- lapply(nam, function(feature)
		{
			elem <- object2[feature][[1]]
			assert.is(elem, "inverseCDF")
			feature.i <- which(nam == feature)
#yo comente 18/10/2012			#if(feature %in% progress.nam)
#yo comente 18/10/2012			#	message("  adding point mass to ", feature, ", ", feature.i, " of ", length(nam), " on ", date())
			addPointMass(elem, mass0 = mass0[feature], param0 = param0, call.browser = call.browser, ...)
		})
		object2@.Data <- lis
		object2@mass0 <- Numeric(mass0[nam])
		stopifnot(validObject(object2))
		if(is.nothing(object2@mass0))
		{ message("bad mass0"); browser()}
		object2
	}
	else
		stop("Cannot add the point mass.")
}

setMethod("p.value", signature(object = "numeric"), function(object, mass0 = NULL, convert.to.two.sided = default(TRUE, "convert.to.two.sided"))
{
	p <- object
	if(is.nothing(mass0))
	{
		mass0 <- numeric(0)
	}
	else if(!is.nothing(mass0) && is(mass0, "empiricalNull"))
	{
		message("getting local false discovery rate estimates on ", date())
		mass0 <- Mass0(mass0)
	}
	else
		stopifnot(!is.nothing(mass0) && is.prob(mass0))
	stopifnot(is(mass0, "numeric") && is.prob(p))
	p2 <- if(convert.to.two.sided)
		2 * pmin(p, 1 - p)
	else
		p
	names(p2) <- names(p)
	stopifnot(is.prob(p2))
	p2
	rm(p)
	if(is.nothing(mass0))  #is.nothing(mass0) && P0 == 0)
	{
		p2
	}
	else if(is.prob(mass0))
	{
		message("invalid p-values (h101225g)"); browser()
		assert.is(mass0, "numeric")
		if(length(mass0) == 1)
		{
			mass0 <- rep(mass0, length(p2))
			names(mass0) <- names(p2)
		}
		stopifnot(is.prob(mass0))
		nam <- intersect(names(mass0), names(p2))
		p2 <- p2[nam]
		mass0 <- mass0[nam]
		mass1 <- 1 - mass0
		stopifnot(length(mass1) == length(p2))
		mass0 * 1 + mass1 * p2
	}
	else
		stop("cannot get p-value")
})
setMethod("p.value", signature(object = "cvalue"), function(object, mass0 = NULL, convert.to.two.sided = default(TRUE, "convert.to.two.sided"), ...)
{
	p <- try(congruity(object = object, ...))
	if(is.err(p))
	{ message("eeee"); browser()}
	assert.is(p, "numeric")
	p.value(object = p, mass0 = mass0, convert.to.two.sided = convert.to.two.sided, ...)
})








when.interval.last.loaded <- date()

# EOF
