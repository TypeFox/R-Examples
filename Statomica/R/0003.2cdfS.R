# cdf.s (formerly fiducial.s) created by David Bickel on 20 April 2009.

new.CDF <- function(object, min_param, max_param, param.name, type)
{
	new("CDF", object, min_param = scalar(min_param), max_param = scalar(max_param), param.name = param.name, type = type)
}
CDF <- function(type, ...)
{
	fun <- if(type == "confidence")
		confidence.CDF
	else if(type == "probability")
		probability.CDF
	else
		stop("bad type of CDF")
	fun(...)
}

blank.CDF <- function(object, param.name = "no parameter")
{
	if(missing(object))
		object <- "no function specified for CDF"
	assert.is(object, "character")
	fun <- function(param)
	{
		stop(object)
	}
	new.CDF(object = fun, min_param = -Inf, max_param = Inf, param.name = param.name, type = "no type")
}

confidence.CDF <- function(object, pvalue.fun, param.name, min_param, max_param, ...) # cf. CD.cdf.binom from logical.r; s3.test.fun removed
{
	if(missing(pvalue.fun) && missing(param.name) && missing(min_param) && missing(max_param))
	{
		message("Due to missing arguments, calling t_test_CDF.")
		t_test_CDF(object = object, ...)
	}
	else
	{
		stopifnot(length(min_param) == 1 && length(max_param) == 1)
		if(missing(object))
		{ message("confidence.CDF missing object"); browser()}
		cdf <- function(param) # integrate needs vector param
		{
#			if(length(param) != 1)
#			{ message("bad param length in confidence.CDF"); browser()}
			confid <- try(ifelse(param < min_param, 0, ifelse(param > max_param, 1, pvalue.fun(object = object, param = param, ...)))) # as of 090509, no longer <= (for integrate)
			if(is(confid, "try-error"))
			{ message("Does pvalue.fun allow vector param?"); print(list(object = object, param = param, ...)); browser()}
			confid
		}
		new.CDF(object = cdf, min_param = min_param, max_param = max_param, param.name = param.name, type = "confidence")
	}
}
rnorm.CDF <- function(object, exponent, sd = 1, n = 1e5, plot = TRUE, call.abs = TRUE, use.confidence = FALSE, ...)
{
	assert.is(object, "numeric")
	stopifnot(length(object) == 1)
	assert.is(exponent, "numeric")
	stopifnot(length(exponent) == 1)
	trans <- function(x){(if(call.abs) abs(x) else x) ^ exponent}
	x0 <- rnorm(mean = 0, sd = sd, n = n)
	get.null.trans.x <- function(param)
	{
		if(length(param) != 1)
		{ message("bad param"); browser()}
		trans(x0 + param) # rnorm(mean = param, sd = sd, n = n)
	}
	pre.fun <- function(object, param)
	{
		observed.trans.x <- trans(object)
		stopifnot(length(observed.trans.x) == 1)
		sapply(param, function(param1){mean(get.null.trans.x(param = param1) > observed.trans.x)})
	}
	pvalue.fun <- if(use.confidence)
	{
		if(call.abs)
			stop("incompatible args")
		else
			function(object, param)
			{
				cdf <- function(q){pre.fun(object = object, param = q)}
				stopifnot(param >= 0)
				cdf(param) - cdf(-param)
			}
	}
	else
		pre.fun
	min_param <- 0
	param.name <- paste("|mean|") #, exponent)
	if(plot)
	{
		delta <- 1e-3
		param <- seq(min_param + delta, 4, delta)
		observation.text <- function(observation){if(call.abs) paste("|", observation, "|", sep = "") else observation}
		plot(x = param, y = pvalue.fun(object = object, param = param), xlab = param.name, ylab = if(use.confidence) "confidence" else "p-value", main = paste("p-value function, x = ", observation.text(signif(object, 1)), "^", exponent, sep = ""), ylim = c(0,1))
	}
	confidence.CDF(object = object, pvalue.fun = pvalue.fun, param.name = param.name, min_param = min_param, max_param = Inf, ...)
}
t_test_CDF <- function(object, ...)
{
	nsample <- if(is(object, "list") && length(object) == 2)
		2
	else if(is(object, "numeric"))
		1
	else
		stop(paste("cannot assign nsample from given object of class", class(object)))
	param.name <- if(nsample == 2)
		"difference of means"
	else if(nsample == 1)
		"mean"
	else
		stop(paste("cannot assign param.name from given object of class", class(object)))
	pvalue.fun <- function(object, param)
	{
		get.pvalue <- function(...)
		{
			t.test(..., mu = param, alternative = "greater")$p.value
		}
		if(nsample == 2)
			get.pvalue(x = object[[1]], y = object[[2]])
		else if(nsample == 1)
			get.pvalue(x = object)
		else
			stop(paste("cannot compute", param.name, "p-value from given object of class", class(object)))
	}
	confidence.CDF(object = object, pvalue.fun = pvalue.fun, param.name = param.name, min_param = -Inf, max_param = Inf, ...)
}
z.test.CDF <- function(object, ...)
{
	nsample <- if(is(object, "list") && length(object) == 2)
		2
	else if(is(object, "numeric"))
		1
	else
		stop(paste("cannot assign nsample from given object of class", class(object)))
	param.name <- if(nsample == 2)
		"difference of means"
	else if(nsample == 1)
		"mean"
	else
		stop(paste("cannot assign param.name from given object of class", class(object)))
	pvalue.fun <- function(object, param)
	{
		get.pvalue <- function(x, y)
		{
			if(missing(y))
			{
				x <- x[is.finite(x)]
				pnorm(mean(x), mean = param, sd = 1 / sqrt(length(x)), lower.tail = FALSE)
			}
			else
				stop("not yet")
		}
		if(nsample == 2)
			get.pvalue(x = object[[1]], y = object[[2]])
		else if(nsample == 1)
			get.pvalue(x = object)
		else
			stop(paste("Cannot compute", param.name, "p-value from given object of class", class(object)))
	}
	confidence.CDF(object = object, pvalue.fun = pvalue.fun, param.name = param.name, min_param = -Inf, max_param = Inf, ...)
}
stat.CDF <- function(object, param.name, statfun, rdata, min_param, max_param, nresample, change.sign = FALSE, ...)
{
	assert.is(object, "numeric")
	assert.is(param.name, "character")
	assert.is(statfun, "function")
	assert.is(rdata, "function")
	assert.is(min_param, "numeric")
	assert.is(max_param, "numeric")
	assert.is(nresample, "numeric")
	assert.is(change.sign, "logical")
	object <- object[is.finite(object)]
	get.data <- function(param)
	{
		rdata(param, n = length(object))
	}
	pvalue.fun <- function(object, param)
	{
		pivotal <- statfun(object) # not necessarily a different object from main argument object
		pivotals <- sapply(1:nresample, function(dummy)
		{
			statfun(get.data(param = param))
		})
		stopifnot(length(pivotal) == 1)
		boo <- if(change.sign)
			pivotals < pivotal
		else
			pivotals > pivotal
		sum(boo) / length(boo)
	}
	confidence.CDF(object = object, pvalue.fun = pvalue.fun, param.name = param.name, min_param = min_param, max_param = max_param, ...)
}
boundedMean.statistical.CDF <- function(object, min_param = default(0, "boundedMean.CDF min_param"), sd = default(1, "boundedMean.CDF sd"), ...) # based on "Is Bayes posterior just quick and dirty confidence?"
{
	if(missing(object))
	{ message("boundedMean.statistical.CDF missing object"); browser()}
	param.name <- "bounded mean"
	pvalue.fun <- function(object, param)
	{
		object <- object[is.finite(object)]
		sample.mean <- mean(object)
		stopifnot(length(sd) == 1 && length(sample.mean) == 1)
		cdf <- function(q)
		{
			stopifnot(length(q) %in% c(1, length(param)))
			pnorm(mean = param, sd = sd / sqrt(length(object)), q = q, lower.tail = FALSE)
		}
		if(any(param >= 0))
			cdf(q = sample.mean) # / cdf(q = param)
		else
			0
	}
	confidence.CDF(object = object, pvalue.fun = pvalue.fun, param.name = param.name, min_param = min_param, max_param = Inf, ...)
}

.delta <- 1e-4

probability.CDF <- function(object, param.name = default("parameter", "param.name"), min_param, max_param, ...)
{
	cdf <- if(is(object, "function"))
		object
	else if(is(object, "numeric")) # has randomly generated parameter values
	{
		if(missing(min_param))
			min_param <- min(object)
		if(missing(max_param))
			max_param <- max(object)
		function(param)
		{
			ecdf(object)(param)
		}
	}
	else
		stop("bad object class")
	new.CDF(object = cdf, min_param = min_param, max_param = max_param, param.name = param.name, type = "probability")
}
t_posterior_CDF <- function(object, ndraw, ...)
{
	nsample <- if(is(object, "list") && length(object) == 2)
		2
	else if(is(object, "numeric"))
		1
	else
		stop(paste("t_posterior_CDF cannot assign nsample from given object of class", class(object)))
	rParam <- function(x){rNormalParameter(x)}
#	get.mean <- function(rp){Ran(rp, name = "mean", ndraw = ndraw)}
#	get.var <- function(rp){Ran(rp, name = "var", ndraw = ndraw)}
	get.cdf <- function(param, param.name)
	{
		probability.CDF(object = param, param.name = param.name, min_param = -Inf, max_param = Inf, ...)
	}
	if(nsample == 1)
	{
		if(missing(ndraw))
		{
			cdf <- t_test_CDF(object = object, ...)
			cdf@type <- "probability"
			cdf
		}
		else # for testing
		{
			rp <- rParam(object)
			param <- Ran(rp, name = "mean", ndraw = ndraw)
			print(stats(param))
			if(mean(param) > 10 * mean(object))
			{ message("large param"); browser()}
			get.cdf(param = param, param.name = "mean")
		}
	}
	else if(nsample == 2)
	{
		assert.is(ndraw, "numeric")
		get.rParam <- function(i){rParam(object[[i]])}
		rp1 <- get.rParam(1)
		rp2 <- get.rParam(2)
		param <- Ran.diff(x = rp1, y = rp2, name = "mean", ndraw = ndraw)
		get.cdf(param = param, param.name = "difference of means")
	}
	else
		stop("bad nsample")
}
Lik.fun.norm <- function(x, sd = default(sd(x), "Lik.fun.norm sd"))
{
	x <- x[is.finite(x)]
	if(length(x) == 0) # no data
	{
		function(param){1}
	}
	else
	{
#		sample.mean <- mean(x)
		stopifnot(length(sd) == 1)
		function(param)
		{
			Lik <- sapply(param, function(mean)
			{
				lik <- dnorm(x = x, mean = mean, sd = sd, log = TRUE) # sapply(x, function(q){dnorm(x = q, mean = param, sd = sd, log = TRUE)})
				stopifnot(length(lik) == length(x))
				exp(sum(lik))
			})
			stopifnot(is.numeric(Lik) && length(Lik) == length(param))
			Lik
		}
	}
}

# boundedMean functions moved to posterior.s on 7 July 2010.

# general code moved to posterior.s on 7 July 2010.

# interval.s code moved from cdf.s on 24 December 2010.

when.cdf.last.loaded <- date()

# EOF