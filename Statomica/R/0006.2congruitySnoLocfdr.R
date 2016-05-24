# congruity.s created by David Bickel on 9 July 2009.
# indicatorRisk changed to probabilisticRisk on 090717 15:19
# plot.risks changed to probabilisticRisks on 090717 15:35

new.marginal.zz <- function(object, pre.marginal){new("marginal.zz", object, pre.marginal = pre.marginal)}
marginal.zz <- function(zz1, zz2)
{
	assert.are(list(zz1, zz2), "numeric")
	combine <- function(x, y){c(as(x, "numeric"), as(y, "numeric"))}
	vec <- combine(zz1, zz2)
	if(is(zz1, "marginal.zz") && is(zz2, "marginal.zz"))
	{
		new.marginal.zz(vec, pre.marginal = combine(zz1@pre.marginal, zz2@pre.marginal))
	}
	else if(!is(zz1, "marginal.zz") && !is(zz2, "marginal.zz"))
		vec
	else
		stop("marginal.zz failed")
}

meanFUN <- function(mean, sd, prob, contamination, ...) # for random.scalar
{
	if(missing(prob))
	{
		prob <- if(missing(contamination))
			1
		else if(length(mean) == 3)
			default(c(contamination / 2, 1 - contamination, contamination / 2), "prob")
		else
			stop("cannot assign prob for meanFUN")
	}
	assert.are(list(mean, sd, prob), "numeric")
	if(length(mean) > 1 && length(sd) == 1)
		sd <- rep(sd, length(mean))
	stopifnot(all(sapply(list(mean, sd), length) == length(prob)))
	stopifnot(sum(prob) == 1)
	function()
	{
		normal.variate <- function(mu, sigma)
		{
			stopifnot(length(mu) == 1 && length(sigma) == 1)
			rnorm(n = 1, mean = mu, sd = sigma, ...)
		}
		ran <- if(length(prob) == 1 && prob == 1)
			normal.variate(mu = mean, sigma = sd)
		else if(length(prob) > 1 && all(prob > 0 & prob < 1))
		{
			i <- rindex(prob = prob)
			normal.variate(mu = mean[i], sigma = sd[i])
		}
		else
			stop("bad meanFUN prob arg")
		ok <- is.numeric(ran) && length(ran) == 1 && is.finite(ran)
		if(!ok)
		{ message("bad ran"); browser()}
		scalar(ran)
	}
}
sdFUN <- function(mean, sd, min, max, FUN, ...) # for random.sd
{
	call.exp <- missing(FUN)
	get.ran <- if(missing(mean) && missing(sd) && missing(min) && missing(max) && !call.exp)
		function(){FUN()}
	else if(missing(mean) && missing(sd))
		function(){runif(n = 1, min = min, max = max)}
	else if(missing(min) && missing(max))
		function(){rnorm(n = 1, mean = mean, sd = sd, ...)}
	else
		stop("sdFUN args incompatible")
	function()
	{
		ran <- get.ran()
		Scalar(if(call.exp) exp(ran) else ran)
	}
}
new.random.scalar <- function(fixed, random)
{
	new("random.scalar", fixed = scalar(fixed), random = random)
}
random.scalar <- function(fixed, sd, contamination, ...) # for mean argument of probabilisticRisk
{
	assert.are(list(fixed, sd), "numeric")
	stopifnot(missing(contamination))
	stopifnot(length(fixed) == 1)
	stopifnot(length(sd) == 1)
	random <- meanFUN(mean = fixed, sd = sd, ...)
	new.random.scalar(fixed = fixed, random = random)
}
random.sd <- function(fixed, sd, ratio, FUN, ...) # for sd argument of probabilisticRisk
{
	assert.is(fixed, "numeric")
	stopifnot(length(fixed) == 1)
	random <- if(missing(sd) && missing(ratio) && !missing(FUN) && is.function(FUN))
	{
		sdFUN(FUN = FUN, ...)
	}
	else if(missing(ratio) && is.numeric(sd))
	{
		stopifnot(length(sd) == 1)
		sdFUN(mean = logb(fixed), sd = sd, ...)
	}
	else if(is.numeric(ratio) && missing(sd))
	{
		if(length(ratio) == 1 && ratio >= 1)
		{
			if(ratio == 1)
				warning("unit ratio will make random element act as fixed")
			exp.min <- fixed / ratio
			exp.max <- fixed * ratio
		}
		else if(length(ratio) == 2 && ratio[1] <= ratio[2])
		{
			if(ratio[1] == ratio[2])
				warning("equal ratio components will make random element act as fixed")
			exp.min <- fixed * ratio[1]
			exp.max <- fixed * ratio[2]
		}
		else
			stop("bad ratio arg")
		stopifnot(exp.min <= exp.max)
		sdFUN(min = log(exp.min), max = log(exp.max), ...)
	}
	else
		stop("arg error in random.sd")
	new.random.scalar(fixed = fixed, random = random)
}
random.df <- function(FUN) # for df argument of probabilisticRisk
{
	assert.is(FUN, "function")
	new.random.scalar(fixed = as.numeric(NA), random = FUN)
}
fixed.value <- function(object)
{
	assert.is(object, "random.scalar")
	object@fixed
}
random.value <- function(object)
{
	assert.is(object, "random.scalar")
	scalar(object@random())
}
extended.scalar <- function(object)
{
	sof <- if(is(object, "random.scalar"))
		object
	else
		scalar(object)
	assert.is(sof, "extended.scalar")
	sof
}
value <- function(object, random = TRUE)
{
	val <- if(is(object, "random.scalar"))
	{
		if(random)
			random.value(object)
		else
			fixed.value(object)
	}
	else if(is.numeric(object))
		scalar(object)
	else
		stop("bad object")
	assert.is(val, "scalar", text = "perhaps mean, sd, mean0, or sd0 is a function that does not return a numeric; you might want to use mean = random.scalar(fixed = 2.5, sd = 1)")
	val
}

new.dataGenerator <- function(object, CDF0, pvalue.fun, zz)
{
	assert.is(CDF0, "function")
	if(!is(CDF0, "CDF"))
		CDF0 <- new.CDF(object = CDF0, min_param = -Inf, max_param = Inf, param.name = "zz", type = "probability")
	new("dataGenerator", object, CDF0 = CDF0, pvalue.fun = pvalue.fun, zz = zz)
}
dataGenerator <- function(object, ngeneration = numeric(0), nobservation = numeric(0), pvalue.fun = default(function(x){t.test(x, alternative = "greater")$p.value}, "pvalue.fun"), verbose = FALSE)
{
	assert.is(object, "function")
	CDF0 <- if(length(ngeneration) == 0)
	{
		zz <- numeric(0)
		function(zz){stop("not a null distribution")}
	}
	else
	{
		zz <- sapply(1:ngeneration, function(i)
		{
			x <- if(is.nothing(nobservation))
				object()
			else
				object(n = nobservation)
			stopifnot(sum(is.finite(x)) >= 2)
			pval <- pvalue.fun(x = x)
			qnorm(pval)
		})
		ECDF(zz, rule = 2)
#		function(v)
#		{
#			cdf <- ECDF(zz, rule = 2) # true conditional null distribution given object as the sampling distribution # fun <- ECDF(zz)
#			p <- cdf(v)
#			if(any(is.na(p)))
#			{ message("any(is.na(p))"); print(summary(p)); browser()}
#			p
#		}
	}
	if(verbose && length(ngeneration) == 1)
		plot(CDF0)
	new.dataGenerator(object = object, CDF0 = CDF0, pvalue.fun = pvalue.fun, zz = zz)
}
rStable <- function(alpha, beta = 0, delta = 0, ...)
{
	function(n) # this n will take nobservation in dataGenerator
	{#XXX|:no visible global function definition for ÔrstableÕ, it exists in a package named "stabledist"
		rStable(n = n, alpha = alpha, beta = beta, delta = delta, ...)#rstable
	}
}
StableGenerator <- function(ngeneration, nobservation = numeric(0), verbose = FALSE, ...)
{
	object <- rStable(...)
	stopifnot(length(nobservation) == 1)
	if(verbose)
		message("calling dataGenerator")
	dataGenerator(object = object, ngeneration = ngeneration, nobservation = nobservation, verbose = verbose)
}

pNorm <- function(Q, Mean, Sd, Df = numeric(0))
{
	if(length(Df) == 0)
		pnorm(q = Q, mean = Mean, sd = Sd)
	else if(length(Df) == 1)
	{
		pt(q = (Q - Mean) / Sd, df = Df)
	}
	else
		stop("Df is of the wrong length")
}
rNorm <- function(N, Mean, Sd, Df = numeric(0), Rfun = NULL, nobservation = numeric(0), uncond.cval)
{
	assert.is(N, "numeric")
	if(!missing(Rfun))
		assert.are(list(Mean, Sd), "numeric")
	direct <- length(nobservation) == 0 && is.null(Rfun)
	if(direct) # sim z-values directly
		stopifnot(length(Df) == 0)
	else if(missing(Rfun)) # sim non-reduced data
	{
		stopifnot(length(Df) == 1 && Df >= 1)
		Rfun <- function(n)
		{
			rt(df = Df, n = n) * Sd + Mean
		}
	}
	assert.is(Rfun, "dataGenerator")
	rfun <- if(!direct && length(nobservation) > 0)
	{
		assert.is(Rfun, "function")
		function()
		{
			ran <- try(Rfun(n = nobservation))
			if(is(ran, "try-error"))
			{ message("is not good"); browser()}
			ran
		}
	}
	else
		Rfun
	zz <- if(missing(uncond.cval) || is.null(uncond.cval))
	{
		if(direct) # sim z-values directly
			rnorm(n = N, mean = Mean, sd = Sd)
		else if(length(nobservation) == 0 || nobservation >= 2) # sim non-reduced data; added 100304: indirect computation of p-values
		{
			p <- sapply(1:N, function(dummy)
			{
				x <- rfun()
				stopifnot(length(nobservation) == 0 || length(x) == nobservation)
				Rfun@pvalue.fun(x) # t.test(x, alternative = "greater")$p.value
			})
			qnorm(p)
		}
		else
			stop("you have nobservation problem")
	}
	else if(is(uncond.cval, "unconditional.cvalue") && direct)
	{
		fun <- as(uncond.cval, "function")
		pre.marginal <- rNorm(N = N, Mean = Mean, Sd = Sd, Df = Df, nobservation = nobservation)
		uncond.zz <- fun(pre.marginal) # converts z-score to z-score that is N(0,1) under the true global null
		new.marginal.zz(object = uncond.zz, pre.marginal = pre.marginal)
	}
	else
		stop("bad uncond.cval combined with direct; cannot compute zz")
	stopifnot(length(zz) == N)
	zz
}


new.Scalar.factor <- function(object){new("Scalar.factor", object)}
Scalar.factor <- function(object)
{
	if(is(object, "Scalar.factor"))
		object
	else
	{
		if(!is(object, "Scalar"))
			object <- Scalar(object)
		new.Scalar.factor(object)
	}
}

new.cvalue <- function(pvalue, zz, s3FUN, arglis)
{
	if(missing(arglis))
	{
		arglis <- list(alternative = "greater")
	}
	if(!missing(pvalue))
	{
		assert.is(pvalue, "Numeric")
		stopifnot(all(pvalue <= 1.01, na.rm = TRUE))
		stopifnot(missing(zz))
		zz <- qnorm(as(pvalue, "numeric"))
		names(zz) <- names(pvalue)
	}
	assert.is(zz, "numeric")
	if(is(zz, "Numeric"))
	{ message("bad zz; should be pvalue?"); browser()}
	new("cvalue", zz = zz, s3FUN = s3FUN, arglis = arglis)
}
cvalue <- function(object, s3FUN, alternative, ttest.arglis, verbose = TRUE, ...) # s3FUN, e.g., t.test, has arguments x, alternative, and possibly y; returns list with "pvalue" in names
{
	if(missing(ttest.arglis))
	{
		if(verbose)
			message("\ncomputing cvalue object on ", date())
		if(missing(alternative))
			alternative <- default("Greater", "\n  cvalue alternative", verbose = verbose)
		alternative <- alt(alternative)
		if(missing(s3FUN) || !is.function(s3FUN))
			stop(paste('s3FUN missing or wrong class; supply argument s3FUN, e.g., t.test,', 'which takes arguments x, alternative, and possibly y; returns list with "pvalue" in names'))
		pvalue <- PValue(object = object, FUN = s3FUN, alternative = alternative, verbose = verbose, ...)
		if(!is.na(pvalue[1]) && all(pvalue == pvalue[1], na.rm = TRUE))
		{ message("all p-values the same"); browser()}
		new.cvalue(pvalue = pvalue, s3FUN = s3FUN, arglis = list(alternative = alternative, ...))
	}
	else if(missing(s3FUN) && missing(alternative))
	{
		assert.is(ttest.arglis, "list")
		t.obj <- do.call("ttest", c(list(object), ttest.arglis))
		as(t.obj, "cvalue")
	}
	else
		stop("cvalue arg err")
}

new.ran.cvalue <- function(cval, param.sign, p0, mean, sd, df, mean0, sd0, df0)
{
	new("ran.cvalue", cval, param.sign = param.sign, p0 = Scalar(p0), mean = scalar(mean), sd = Scalar(sd), df = Scalar(df), mean0 = scalar(mean0), sd0 = Scalar(sd0), df0 = Scalar(df0))
}
ran.cvalue <- function(n, n0, nobservation = numeric(0), mean = default(2.5, "mean"), sd = default(1.25, "sd"), df = numeric(0), rfun = NULL, mean0 = default(0, "mean"), sd0 = default(1, "sd0"), df0 = numeric(0), rfun0 = NULL, uncond.cval = NULL, verbose = FALSE) # defaults from Fig. 5 of Efron, B. (2007) "Correlation and Large-Scale Simultaneous Significance Testing", Jour Amer Stat Assoc, 102, pp. 93-103
{
	stopifnot((is(rfun, "function") && is(rfun0, "function")) || (is.null(rfun) && is.null(rfun0)))
	stopifnot(is.null(uncond.cval) || is(uncond.cval, "unconditional.cvalue")) # added 090712
	stopifnot(length(n) == 1)
	direct <- length(nobservation) == 0
	if(direct) # sim z-values directly
		stopifnot(length(df) == 0 && length(df0) == 0)
	else # sim non-reduced data
		stopifnot(length(df) >= 1 || is.function(rfun)) # length(df) == 1 && length(df0) == 1 && df >= 1 && df0 >= 1)
	n1 <- n - n0
	fixed.mean <- value(mean, random = FALSE) # for param.sign
	fixed.mean0 <- value(mean0, random = FALSE) # for param.sign
	mean0 <- value(mean0)
	sd0 <- value(sd0) # error corrected 090711 08:11
	mean <- if(is(mean, "Scalar.factor"))
		mean * sd0 # mean rescaled by possibly random sd0
	else
		value(mean)
	sd <- if(is(sd, "Scalar.factor"))
		sd * sd0
	else
		value(sd)
	df <- value(df)
	df.ok <- function(dof){length(dof) == 0 || (is.finite(dof) && dof >= 1)}
	if(!df.ok(df))
	{ message("bad degree of freedom"); browser()}
	df0 <- if(length(df0) == 0)
		df
	else
		value(df0)
	if(!df.ok(df0))
	{ message("bad degree of freedom 0"); browser()}
	get.file <- function(...){paste(..., " sd0=", if(floor(sd0) == sd0) sd0 else "fraction", ".txt", sep = "")}
	print2file("cvalue:  mean", as(mean, "numeric"), "sd", as(sd, "numeric"), file = get.file("cvalue mean and sd"))
	get.zz <- function(N, Mean, Sd, Df, Rfun)
	{
		rNorm(N = N, Mean = Mean, Sd = Sd, Df = Df, Rfun = Rfun, nobservation = nobservation, uncond.cval = uncond.cval)
	}
	zz1 <- get.zz(N = n1, Mean = mean, Sd = sd, Df = df, Rfun = rfun)
	zz2 <- get.zz(N = n0, Mean = mean0, Sd = sd0, Df = df0, Rfun = rfun0)
	zz <- marginal.zz(zz1, zz2) # order must match coverage function and param.sign
	print2file(stats(zz[1:n1]), FUN = print, file = get.file("cvalue first n1 of zz"))
	print2file(stats(zz[(n1 + 1):length(zz)]), FUN = print, file = get.file("cvalue last n0 of zz"))
	get.param.sign <- function(nsign, z.mean)
	{
		if(z.mean == 0)
			sign(runif(nsign) - 0.5) # there is no information in c-value about the sign
		else
			rep(-sign(z.mean), nsign) # a high c-value means a high congruity of the data with a parameter smaller than the null value
	}
	param.sign <- c(get.param.sign(nsign = n1, z.mean = fixed.mean), get.param.sign(nsign = n0, z.mean = fixed.mean0)) # added 090710
	stopifnot(length(param.sign) == length(zz))
	if(length(zz) != n)
	{ message("bad zz"); browser()}
	if(xor(!is.null(uncond.cval), is(zz, "marginal.zz")))
	{ message("zz incompatible with uncond.cval"); browser()}
	cval <- new.cvalue(zz = zz, s3FUN = ran.cvalue, arglis = list(alternative = "Greater", n = n, n0 = n0))
	new.ran.cvalue(cval = cval, param.sign = param.sign, p0 = n0 / n, mean = mean, sd = sd, df = df, mean0 = mean0, sd0 = sd0, df0 = df0)
}

new.adjusted.cvalue <- function(object, empirical.null)
{
	new("adjusted.cvalue", object, empirical.null = empirical.null)
}
adjusted.cvalue <- function(object, empirical.null, CDF0, adjust.pre.marginal, allow.class.change = TRUE, ...) # CDF0 argument added 090717 for conditional.con of probabilisticRisk; adjust.pre.marginal added 090721 to correct error with marginally correct c-values
{
	if(missing(object))
	{ message("object arg missing"); browser()}
	assert.is(object, "cvalue")
	if(!xor(missing(empirical.null), missing(CDF0)))
	{ message("exactly one of these args should be specified: empirical.null; CDF0"); browser()}
	has.empirical.null <- !missing(empirical.null)
	if(missing(empirical.null) && missing(CDF0))
	{
		empirical.null <- empiricalNull(object = object, ...)
		has.empirical.null <- TRUE
	}
	if(has.empirical.null)
		assert.is(empirical.null, "empiricalNull")
	nominal.zz <- if(has.empirical.null && is.character(names(empirical.null)))
		object@zz[names(empirical.null)]
	else
		object@zz
	if(missing(CDF0))
		CDF0 <- empirical.null@CDF0
	assert.is(CDF0, "CDF")
	if(missing(adjust.pre.marginal))
		adjust.pre.marginal <- !has.empirical.null && is(nominal.zz, "marginal.zz")
	p <- sapply(if(adjust.pre.marginal) nominal.zz@pre.marginal else nominal.zz, CDF0) # h100429
	if(!is.numeric(p))
	{ message("bad p of adjusted.cvalue"); browser()}
	object@zz <- qnorm(p = p) # [names(empirical.null)] added 090716 to pass validObject
	acval <- if(has.empirical.null)
		new.adjusted.cvalue(object = object, empirical.null = empirical.null)
	else
		try(object)
	if(!allow.class.change && class(acval) != class(object))
	{
		acval <- if(is(object, "statistic.cvalue"))
		{
			object@zz <- acval@zz
			object@statistic <- object@statistic[names(object@zz)]
			object
		}
		else if(is(object, "cvalue"))
			as(acval, "cvalue")
		else
			stop("class could not be retained")
	}
	if(!is(acval, "cvalue") || !validObject(acval))
	{ message("adjustment error"); browser()}
	acval
}

new.statistic.cvalue <- function(object, statistic, statistic.name, statistic.fun)
{
	new("statistic.cvalue", object, statistic = statistic, statistic.name = statistic.name, statistic.fun = statistic.fun)
}
statistic.cvalue <- function(object, statistic.name, statistic.fun, ...)
{
	if(missing(statistic.name) && missing(statistic.fun))
	{
		statistic.name <- default("mean", "statistic.name")
		statistic.fun <- default(mean, "statistic.fun")
	}
	assert.is(statistic.name, "character")
	assert.is(statistic.fun, "function")
	statistic <- stat(object, FUN = statistic.fun)
	cval <- cvalue(object = object, ...)
	if(length(statistic) != length(cval))
	{ message("bad statistic"); browser()}
	new.statistic.cvalue(object = cval, statistic = statistic, statistic.name = statistic.name, statistic.fun = statistic.fun)
}

default.sign.loss.fun <- function(true.sign, decision, npass, b = default(1.5, "default.sign.loss.fun b"), a = 1, c = 0)
{
	stopifnot(length(true.sign) == length(decision))
	stopifnot(all(true.sign %in% c(-1, 1)))
	stopifnot(all(decision %in% c(-1, 0, 1)))
	if(missing(npass))
		npass <- sum(decision == 0, na.rm = TRUE) # was true.sign
	nwrong <- sum(decision != 0 & decision != true.sign)
	loss.term <- function(nwrong)
	{
		if(length(b) == 0 || b < 0)
			a * (nwrong + c)
		else if(length(b) == 1)
			a * (nwrong + c) ^ b # a * b ^ (nwrong + c)
		else
			stop("bad b")
	}
	guessing.loss <- loss.term(nwrong = nwrong) - loss.term(nwrong = 0)
	guessing.loss + npass
}
sign_expected_loss <- function(prob.negative, decision, sign.loss.fun, nsample = default(100, "sign_expected_loss nsample"), b = default(1.5, "sign_expected_loss b", verbose = missing(sign.loss.fun)), a = 1, c = 0, verbose = FALSE)
{
	if(missing(sign.loss.fun))
		sign.loss.fun <- function(true.sign, decision, npass){default.sign.loss.fun(true.sign = true.sign, decision = decision, npass = npass, a = a, b = b, c = c)}
	assert.are(list(prob.negative, decision), "numeric")
	assert.is(sign.loss.fun, "function")
	if(length(decision) != length(prob.negative))
		stop("mismatch 1")
	boo <- is.finite(decision) & is.finite(prob.negative)
	decision <- decision[boo]
	prob.negative <- prob.negative[boo]
	if(length(decision) != length(prob.negative))
		stop("mismatch 2")
	stopifnot(all(decision %in% (c(-1, 0, +1))))
	decided <- decision != 0
	npass <- sum(!decided)
	prob.negative <- prob.negative[decided]
	decision <- decision[decided]
	if(length(decision) != length(prob.negative))
		stop("mismatch 3")
	if(!all(decision %in% c(-1, +1)))
	{ message("Bad decision!"); browser()}
#	expected.nwrong <- sum(ifelse(decision == -1, 1 - prob.negative, prob.negative))
#	len.ok <- length(decision) >= expected.nwrong && expected.nwrong >= 0
	sel <- if(length(b) == 1 && b >= 0)
	{
		if(verbose)
			message("  conducting ", nsample, " Monte Carlo iterations")
		get.true.sign <- function()
		{
			truth <- ifelse(rbinom(n = length(prob.negative), size = 1, prob = prob.negative), -1, +1)
			if(!all(truth %in% c(-1, +1)))
			{ message("Bad truth!"); browser()}
			truth
		}
		losses <- sapply(1:nsample, function(i){sign.loss.fun(true.sign = get.true.sign(), decision = decision, npass = npass)})
		if(all(losses == losses[1]) && losses[1] != npass)
		{ message("all losses the same"); browser()}
		mean(losses)
	}
	else if(length(b) == 0 || b < 0)
	{
		if(verbose)
			message("  bypassed Monte Carlo; ignoring c")
		eloss <- function(prob.wrong, loss.if.wrong = a)
		{
			if(any(prob.wrong[is.finite(prob.wrong)] > 1/2))
			{ message("foolish decision"); browser()}
			prob.wrong * loss.if.wrong
		}
		elosses <- if(length(decision) > 0)
		{
			prob.wrong <- ifelse(decision == -1, 1 - prob.negative, prob.negative) # ifelse(decision == 0, NA, )
			eloss(prob.wrong = prob.wrong) # ifelse(decision == 0, 1, eloss(prob.wrong = prob.wrong))
		}
		else
			0
		stopifnot(all(is.finite(elosses)))
		sum(elosses) + npass
	}
	else
		stop("bad b length")
	if(any(is.nan(sel)))
	{ message("bad sel"); browser()}
	stopifnot(length(sel) == 1)
	sel
}
sign_decision <- function(object, nsample = default(100, "sign_decision nsample"), save.time = TRUE, b = default(1.5, "sign_decision b", verbose = TRUE), a = default(1, "sign_decision a", verbose = TRUE), c = default(0, "sign_decision c", verbose = TRUE), verbose = FALSE, ...)
{
	assert.is(object, "statistic.cvalue")
	prob.less <- congruity(object)
	decisiveness <- abs(prob.less - 0.5)
	decision0 <- rep(0, length(prob.less))
	nam <- names(prob.less)
	stopifnot(is(nam, "character"))
	stat.sign <- -sign(prob.less - 0.5) # sign(object@statistic[nam])
	make.decision <- function(decisiveness.threshold)
	{
		oks <- c(length(decisiveness) == length(stat.sign), length(decisiveness.threshold) == 1, length(decision0) == length(stat.sign))
		if(!all(oks==T))
		{ print(oks); message("bad decision"); browser()}
		decision <- ifelse(decisiveness >= decisiveness.threshold, stat.sign, decision0)
		stopifnot(length(decision) == length(nam))
		names(decision) <- nam
		if(length(b) == 0 || b < 0)
		{
			prob.wrong <- ifelse(decision == 0, NA, ifelse(decision == -1, 1 - prob.less, prob.less))
			if(any(prob.wrong > 1/2, na.rm = TRUE))
			{ message("Foolishness."); browser()}
		}
		decision
	}
	expected.loss <- function(decisiveness.threshold)
	{
		sign_expected_loss(prob.negative = prob.less, decision = make.decision(decisiveness.threshold), nsample = nsample, a = a, b = b, c = c, verbose = verbose, ...) # sign.loss.fun = sign.loss.fun,
	}
	udecisiveness <- unique(decisiveness[is.finite(decisiveness)])
	decisiveness.threshold <- if(save.time)
	{
#		stop("time-efficient version not yet implemented")
		thresholds <- sort(udecisiveness, decreasing = TRUE)
		candidate.eloss <- eloss <- candidate.threshold <- Inf
		i <- 0
		while(candidate.eloss <= eloss)
		{
			eloss <- candidate.eloss
			threshold <- candidate.threshold
			i <- i + 1
			if(verbose) message("  iteration ", i, " of up to ", length(thresholds), "; threshold: ", signif(threshold, 5), "; eloss: ", signif(eloss, 5))
			if(i > length(thresholds))
				break()
			candidate.threshold <- thresholds[i]
			candidate.eloss <- expected.loss(decisiveness.threshold = candidate.threshold)
		}
		if(threshold < 0)
		{ message("bad threshold"); browser()}
		threshold
	}
	else
	{
		expected.losses <- sapply(udecisiveness, function(decisiveness.threshold)
		{
			expected.loss(decisiveness.threshold = decisiveness.threshold)
		})
		opt.decisiveness <- udecisiveness[expected.losses == min(expected.losses, na.rm = TRUE)]
		min(opt.decisiveness)
	}
	message("Expected loss ", signif(expected.loss(decisiveness.threshold = decisiveness.threshold), 5), " achieved by threshold ", signif(decisiveness.threshold, 3), " on ", date(), ".")
	object@statistic <- make.decision(decisiveness.threshold = decisiveness.threshold)
	object@zz <- object@zz[names(object@statistic)]
	object@statistic.fun <- sign_decision
	object@statistic.name <- "sign decision"
	if(!validObject(object))
	{ message("sign_decision error"); browser()}
	object
}

new.cvalues <- function(object, b, b.name){new("cvalues", object, b = b, b.name = b.name)}
cvalues <- function(object, b, b.name = "risk aversion", nsample = default(1e3, "nsample"), ...)
{
	if(is(object, "list"))
	{
		assert.is(b, "numeric")
		assert.is(b.name, "character")
		len.ok <- length(b) %in% c(0, length(object)) && length(b.name) == 1
		if(!len.ok)
		{ message("bad len in cvalues"); browser()}
		new.cvalues(object = object, b = b, b.name = b.name)
	}
	else if(is(object, "statistic.cvalue"))
	{
		get.sign_decision <- function(B){sign_decision(object, b = B, nsample = nsample, ...)}
		lis <- if(length(b) == 0) # || b < 0)
		{
			B <- numeric(0)
			b <- 1
			list(get.sign_decision(B = B))
		}
		else
			lapply(b, function(B){message("\ndeciding ", B, " on ", date()); get.sign_decision(B = B)})
		cvalues(lis, b = b, b.name = b.name)
	}
	else
		stop("bad args")
}

new.coverage <- function(ncovered, n0, nominal.rate, p0, mean, sd, mean0, sd0, nulltype)
{
	new("coverage", ncovered, n0 = Scalar(n0), nominal.rate = Scalar(nominal.rate), p0 = Scalar(p0), mean = scalar(mean), sd = Scalar(sd), mean0 = scalar(mean0), sd0 = Scalar(sd0), nulltype = Scalar(nulltype))
}
coverage <- function(nsample, n = default(3000, "n"), p0 = 0.95, n0, mean = default(2.5, "mean"), sd = default(1.25, "sd"), mean0 = default(0, "mean0"), sd0 = default(1, "sd0"), nulltype = default(1, "nulltype"), plot = default(0, "plot"), nominal.rate = default(0.95, "nominal.rate"), ...) # defaults from Fig. 5 of Efron, B. (2007) "Correlation and Large-Scale Simultaneous Significance Testing", Jour Amer Stat Assoc, 102, pp. 93-103
{
	stopifnot(length(n) == 1)
	stopifnot(length(nominal.rate) == 1)
	if(missing(n0))
		n0 <- default(round(p0 * n), "n0")
	message("\n\nbeginning simulations on ", date())
	ncovered <- sapply(1:nsample, function(i)
	{
		cval <- ran.cvalue(n = n, n0 = n0, mean = mean, sd = sd, mean0 = mean0, sd0 = sd0)
		acval <- adjusted.cvalue(cval, nulltype = nulltype, plot = plot, ...)
		p <- congruity(acval)[(n - n0 + 1):length(acval)] # order must match ran.cvalue function
		if(length(p) != n0)
		{ message("bad p length"); browser()}
		half.alpha <- (1 - nominal.rate) / 2
		stopifnot(length(half.alpha) == 1)
		covered <- p <= 1 - half.alpha & p > half.alpha
		sum(covered)
	})
	message("ended simulations on ", date())
	stopifnot(length(ncovered) == nsample)
	new.coverage(ncovered = ncovered, n0 = n0, nominal.rate = nominal.rate, p0 = n0 / n, mean = mean, sd = sd, mean0 = mean0, sd0 = sd0, nulltype = nulltype)
}

new.cvalueError.FUN <- function(object, ann, min_congruity_accept, max_congruity_accept)
{
	new("cvalueError.FUN", object, ann = ann, min_congruity_accept = Scalar(min_congruity_accept), max_congruity_accept = Scalar(max_congruity_accept))
}
conservatism <- function(min_congruity_accept, max_congruity_accept) # conservatism of nominal c-values over conditional c-values, i.e., extent to which nominal c-values accept null hypotheses rejected by conditional c-values minus the reverse
{
	object <- function(nominal.congruity, conditional.congruity)
	{
		accept.null <- function(con)
		{
			con >= min_congruity_accept & con <= max_congruity_accept
		}
		accepted.nominally <- accept.null(nominal.congruity)
		accepted.conditionally <- accept.null(conditional.congruity)
		ntests <- length(nominal.congruity)
		lens <- sapply(list(nominal.congruity, conditional.congruity, accepted.nominally, accepted.conditionally), length)
		if(!all(ntests == lens))
		{ message("bad length"); print(lens); browser()}
		nconservative <- sum(accepted.nominally & !accepted.conditionally) - sum(accepted.conditionally & !accepted.nominally) # sign error corrected 090720 12:01; see "pre conservatism sign error corrected.zip"
		stopifnot(length(nconservative) == 1)
		normalized <- nconservative / ntests
		stopifnot(all(abs(normalized) <= 1.0001))
		normalized
	}
	new.cvalueError.FUN(object = object, ann = "conservatism", min_congruity_accept = min_congruity_accept, max_congruity_accept = max_congruity_accept)
}
cvalueError.FUN.slot <- function(object, name)
{
	if(is(object, "cvalueError.FUN"))
		slot(object, name = name)
	else if(is(object, "probabilisticRisk") && is(object@loss.FUN, "cvalueError.FUN"))
		cvalueError.FUN.slot(object = object@loss.FUN, name = name)
	else if(is(object, "probabilisticRisks"))
	{
		vec <- sapply(object, cvalueError.FUN.slot, name = name)
		vec1 <- vec[1]
		stopifnot(all(vec1 == vec))
		vec1
	}
	else
	{ message("cannot extract slot ", name); browser()}
}
min_congruity_accept <- function(object){cvalueError.FUN.slot(object = object, name = "min_congruity_accept")}
max_congruity_accept <- function(object){cvalueError.FUN.slot(object = object, name = "max_congruity_accept")}
setMethod("annotation", signature(object = "probabilisticRisk"), function(object)
{
	cvalueError.FUN.slot(object = object, name = "ann")
})
setMethod("annotation", signature(object = "probabilisticRisks"), function(object)
{
	cvalueError.FUN.slot(object = object, name = "ann")
})
setMethod("export", signature(object = "probabilisticRisks"), function(object, ...)
{
	datf <- do.call("rbind", lapply(object, as, Class = "data.frame"))
	invisible(datf)
})
export.probabilisticRisksLis <- function(object, file, files, ...)
{
	assert.is(object, "list")
	assert.are(object, "probabilisticRisks")
	assert.is(file, "character")
	assert.is(files, "character")
	stopifnot(length(object) == length(files))
	datfs <- lapply(1:length(object), function(i)
	{
		da <- plot(object[[i]], file = files[i], ...)
#		da$file <- files[i]
		nam <- if(length(files[i]) == 1)
			rep(files[i], nrow(da))
		else
			files[i]
		nam <- make.names(nam, unique = TRUE)
		if(nrow(da) != length(nam))
		{ message("nam err"); browser()}
		filenames <- try(rownames(da) <- nam)
		if(is(filenames, "try-error"))
		{ message("filenames err"); browser()}
		da
	})
	datf <- do.call("rbind", datfs)
	write.csv(datf, file = file)
	datf
}

new.probabilisticRisk <- function(cvalue.risk, fdr.risk, n0, p0, mean, sd, mean0, sd0, nulltype, loss.FUN)
{
	new("probabilisticRisk", cvalue.risk = cvalue.risk, fdr.risk = fdr.risk, n0 = Scalar(n0), p0 = Scalar(p0), mean = extended.scalar(mean), sd = extended.scalar(sd), mean0 = extended.scalar(mean0), sd0 = extended.scalar(sd0), nulltype = Scalar(nulltype), loss.FUN = loss.FUN)
}
probabilisticRisk <- function(nsample, n = default(3000, "n"), nobservation = numeric(0), ngeneration = numeric(0), p0 = 0.95, n0, mean = default(2.5, "mean"), sd = default(1.25, "sd"), df = numeric(0), alpha = numeric(0), delta = numeric(0), rfun = NULL, mean0 = default(0, "mean0"), sd0 = default(1, "sd0"), df0 = numeric(0), alpha0 = numeric(0), delta0 = numeric(0), rfun0 = NULL, nulltype = default(1, "nulltype"), plot = default(0, "plot"), loss.FUN, uncond.cval = NULL, verbose = FALSE, congruity.range.accept, max.abs.cvalue.error1 = Inf, nulltype.monitored = c(0, 1), sample.monitored = -1, always.error = FALSE, ...) # defaults from Fig. 5 of Efron, B. (2007) "Correlation and Large-Scale Simultaneous Significance Testing", Jour Amer Stat Assoc, 102, pp. 93-103
{
	get.rfun <- function(a, d, ngeneration)
	{
		if(!is.nothing(a) && !is.nothing(d))
		{
			assert.are(list(a, d), "numeric")
			StableGenerator(ngeneration = ngeneration, nobservation = nobservation, alpha = a, beta = 0, delta = d, verbose = verbose)
		}
		else
		{
			stopifnot(is.nothing(a) && is.nothing(d))
			NULL
		}
	}
	if(is.null(rfun))
		rfun <- get.rfun(a = alpha, d = delta, ngeneration = numeric(0))
	else
	{
		stopifnot(is.nothing(alpha) && is.nothing(delta))
		assert.is(rfun, "function")
	}
	if(is.null(rfun0))
		rfun0 <- get.rfun(a = alpha0, d = delta0, ngeneration = ngeneration)
	else
	{
		stopifnot(is.nothing(alpha0) && is.nothing(delta0))
		assert.is(rfun0, "dataGenerator")
	}
	stopifnot((is(rfun, "function") && is(rfun0, "dataGenerator")) || (is.null(rfun) && is.null(rfun0)))
	if(is(rfun0, "dataGenerator") && !is.nothing(alpha0) && !is.nothing(delta0))
	{
		main <- paste("alpha0=", alpha0, " delta0=", delta0, sep = "") # "dataGenerator CDF0"
		pdf(file = paste(paste(main, "--", n, " features", sep = ""), "pdf", sep = "."))
		plot(rfun0, main = main, congruity.range.accept = congruity.range.accept)
		dev.off()
	}
	stopifnot(length(n) == 1)
	if(missing(loss.FUN))
	{
		loss.FUN <- if(missing(congruity.range.accept) || is.null(congruity.range.accept))
			squaredError
		else if(is.numeric(congruity.range.accept) && length(congruity.range.accept) == 2 && congruity.range.accept[1] <= congruity.range.accept[2])
			conservatism(min_congruity_accept = congruity.range.accept[1], max_congruity_accept = congruity.range.accept[2])
		else
			stop("bad args to probabilisticRisk")
		default(loss.FUN, "loss.FUN")
	}
	assert.is(loss.FUN, "function")
	stopifnot(is.null(uncond.cval) || is(uncond.cval, "unconditional.cvalue"))
	if(missing(n0))
		n0 <- default(round(p0 * n), "n0")
	cvalue.name <- "cvalue"
	fdr.name <- "fdr"
	message.i <- round(seq(1, nsample, length.out = 10))
	message("\n\nbeginning simulations on ", date(), ". message.i:")
	print(message.i)
	risk.lis <- lapply(1:nsample, function(i)
	{
		success <- FALSE
		while(!success)
		{
			cval <- ran.cvalue(n = n, n0 = n0, mean = mean, nobservation = nobservation, sd = sd, df = df, rfun = rfun, mean0 = mean0, sd0 = sd0, df0 = df0, rfun0 = rfun0, verbose = verbose, uncond.cval = uncond.cval)
			assert.is(cval, "ran.cvalue")
			alternative <- alt(cval)
			stopifnot(length(alternative) == 0 || alternative == "Greater")
			en <- empiricalNull(object = cval, nulltype = nulltype, plot = plot, ...)
			success <- is(en, "empiricalNull") && length(cval) == length(en)
			if(!success)
				message("  Sample failed to generate acceptable empirical null; retrying ", date())
		}
		theoretical.null <- cval@mean0 == 0 && cval@sd0 == 1
		if(!is(en, "empiricalNull"))
		{
			return(list())
		}
		fdr <- lfdr(en)
		get.con <- function(...)
		{
			con <- congruity(adjusted.cvalue(cval, ...))
			if(length(con) != n)
			{ message("con err"); browser()}
			con
		}
		acon <- get.con(empirical.null = en) # cval adjusted for estimated null, if any
		get.CDF0 <- function(empirical.null)
		{
			CDF0 <- empirical.null@CDF0
			CDF0@.Data <- function(q){pNorm(Q = q, Mean = cval@mean0, Sd = cval@sd0, Df = cval@df0)} # zz cumulative null distribution conditional on mean0, sd0, df0 that generated the data
			stopifnot(validObject(CDF0) && !identical(CDF0, empirical.null@CDF0))
			CDF0
		}
		if(xor(!is.null(uncond.cval), is(cval@zz, "marginal.zz")))
		{ message("cval@zz incompatible with uncond.cval"); browser()}
		CDF0 <- if(is.null(rfun0))
			get.CDF0(empirical.null = en) # this CDF0 is true conditional null distribution given mean0, sd0, df0
		else if(is(rfun0, "dataGenerator"))
		{
			stopifnot(length(ngeneration) == 1 && ngeneration > 1) # otherwise CDF0 not available
			rfun0@CDF0 # this CDF0 is true conditional null distribution given rfun0 as the sampling distribution
		}
		else
			stop("rfun0 err")
		conditional.con <- get.con(CDF0 = CDF0)
		if(length(conditional.con) != n)
		{ message("conditional.con err"); browser()}
		ok <- !always.error || !identical(conditional.con, acon)
		if(!ok)
		{ message("identical(conditional.con, acon)"); browser()}
		if(length(fdr) != length(acon))
		{
			message("fdr has ", length(fdr), " but acon has ", length(acon))
			browser()
		}
		extract.vec <- function(whole.vec, null)
		{
			if(length(whole.vec) != n)
			{ message("whole.vec err"); browser()}
			n1 <- n - n0
			alt.i <- 1:n1
			null.i <- (n1 + 1):length(whole.vec)
			stopifnot(all(c(alt.i, null.i) == 1:n))
			p <- whole.vec[if(null) null.i else alt.i] # order must match ran.cvalue function
			if((null && length(p) != n0) || (!null && length(p) != n1))
			{ message("bad p length in extract.vec"); browser()}
			p
		}
		get.error <- function(whole.fdr, null)
		{
			whole.prob.neg <- acon
			if(!missing(whole.fdr))
				whole.prob.neg <- ifelse(whole.prob.neg > 0.5, 1 - whole.fdr / 2, whole.fdr / 2) # converts fdr to P(negative) per h090710
				# redundant: ifelse(whole.fdr == 1, 0.5, )
			extract <- function(whole.vec){extract.vec(whole.vec = whole.vec, null = null)}
			if(length(whole.prob.neg) != n)
			{ message("whole.prob.neg err"); browser()}
			prob.neg <- extract(whole.prob.neg)
			err <- if(is(loss.FUN, "cvalueError.FUN"))
			{
				conditional.prob.neg <- extract(conditional.con)
				if(length(prob.neg) != length(conditional.prob.neg))
				{ message("bad conditional.prob.neg length"); browser()}
				prob.ok <- !always.error || !identical(all.equal(conditional.prob.neg, prob.neg), TRUE) # compares conditional.con to acon
				if(!prob.ok)
				{ message("all.equal(conditional.prob.neg, prob.neg)"); browser()}
				loss.FUN(nominal.congruity = prob.neg, conditional.congruity = conditional.prob.neg)
			}
			else
			{
				if(length(cval@param.sign) != n)
				{ message("cval@param.sign err"); browser()}
				neg.indicator <- extract(ifelse(cval@param.sign < 0, 1, 0))
				if(length(prob.neg) != length(neg.indicator))
				{ message("bad neg.indicator length"); browser()}
				loss.FUN(predicted = prob.neg, observed = neg.indicator)
			}
			if(always.error && all(err == 0))
			{ message ("Is there really no err? Come on."); browser()}
			err
		}
		cvalue.error1 <- get.error(null = FALSE)
		if(nulltype %in% nulltype.monitored && (i %in% sample.monitored || abs(cvalue.error1) > max.abs.cvalue.error1))
		{
			sub <- paste(class(uncond.cval), pwd())
			plot.title <- function(){title(main = paste("cvalue.error1:", signif(cvalue.error1, 2)))}
			plot.hist <- function(...){hist(..., breaks = 10); title(sub = sub)}
			get.xlab <- function(mu, sigma){paste("z adjusted by N(", mu, ",", signif(sigma, 2), "^2)", sep = "")}
			lim <- c(-3, 3)
			plot.signif <- function(lim = c(-5, 5))
			{
				plot(qnorm(acon), qnorm(conditional.con), xlab = acon.lab, ylab = ccon.lab, xlim = lim, ylim = lim)
				abline(0, 1, col = "blue")
				plot.title(); title(sub = sub)
				z.acc <- qnorm(congruity.range.accept)
				abline(v = z.acc, h = z.acc, col = "orange")
			}
			par(mfrow = c(2,2))
			en.test <- empiricalNull(object = cval, nulltype = nulltype, plot = 1, ...); plot.title()
			if(nulltype == 0)	plot.signif()	else plot.hist(cval@zz, xlim = lim, xlab = "z from rNorm")
			acon.lab <- get.xlab(mu = Mean(en, nulltype = nulltype), sigma = Sd(en, nulltype = nulltype))
			plot.hist(qnorm(acon), xlim = lim, xlab = acon.lab)
			ccon.lab <- get.xlab(mu = cval@mean0, sigma = cval@sd0)
			plot.hist(qnorm(conditional.con), xlim = lim, xlab = ccon.lab)
			message("cvalue.error1 ", if(is(loss.FUN, "cvalueError.FUN")) paste("(", annotation(loss.FUN), ") ", sep = "") else "", signif(cvalue.error1, 3), " for nulltype ", nulltype, " at sample ", i)
			if(nulltype > 0)
			{
				nx11()
				par(mfrow = c(2,2))
				plot.signif()
#				plot(cval@zz, qnorm(acon), xlim = lim, ylim = lim); plot.title(); title(sub = sub)
				hist(cval@zz@pre.marginal, breaks = 200); title(sub = sub)
				hist(cval@zz, breaks = 200); title(sub = sub)
				hist(qnorm(acon), breaks = 200); title(sub = sub)
			}
			browser()
		}
		fdr.error1 <- get.error(whole.fdr = fdr, null = FALSE)
		cvalue.error0 <- get.error(null = TRUE)
		if(always.error && all(cvalue.error0 == 0))
		{ message ("Is there really no error?"); browser()}
		fdr.error0 <- get.error(whole.fdr = fdr, null = TRUE)
		get.riskPair <- function(error1, error0)
		{
			riskPair(error1 = error1, error0 = error0, loss.FUN = loss.FUN)
		}
		if(i %in% message.i)
			message("    ending iteration ", i, " on ", date())
		lis <- list(get.riskPair(error1 = cvalue.error1, error0 = cvalue.error0), get.riskPair(error1 = fdr.error1, error0 = fdr.error0)) # risks for one trial
		names(lis) <- c(cvalue.name, fdr.name)
		lis
	}) # end lapply i
	message("ended simulations on ", date())
	stopifnot(length(risk.lis) == nsample)
	risk <- function(prob.name)
	{
		stopifnot(prob.name %in% c(cvalue.name, fdr.name))
		get.risk <- function(name)
		{
			vec <- sapply(risk.lis, function(lis)
			{
				assert.is(lis, "list")
				if(length(lis) == 0) # locfdr failed
					as.numeric(NA)
				else
				{
					pair <- lis[prob.name][[1]]
					assert.is(pair, "riskPair")
					slot(pair, name = name)
				}
			})
			if(length(vec) != nsample)
			{ message("bad vec length"); browser()}
			vec
		}
		risk1 <- get.risk("risk1")
		risk0 <- get.risk("risk0")
		new.riskPair(risk1 = risk1, risk0 = risk0, loss.FUN = loss.FUN)
	}
	cvalue.risk <- risk(cvalue.name)
	fdr.risk <- risk(fdr.name)
#	stopifnot(length(cvalue.risk) == nsample)
#	stopifnot(length(fdr.risk) == nsample)
	new.probabilisticRisk(cvalue.risk = cvalue.risk, fdr.risk = fdr.risk, n0 = n0, p0 = n0 / n, mean = mean, sd = sd, mean0 = mean0, sd0 = sd0, nulltype = nulltype, loss.FUN = loss.FUN)
}

default.riskPairFUN <- function(x)
{
	na.boo <- is.na(x)
	if(all(na.boo))
	{ message("no risks to average"); browser()}
	else if(any(na.boo))
		warning(paste(sum(na.boo), "of", length(x), "missing values removed to average risk"))
	mean(x, na.rm = TRUE)
}
new.riskPair <- function(risk1, risk0, FUN = default.riskPairFUN, loss.FUN)
{
	new("riskPair", risk1 = risk1, risk0 = risk0, FUN = FUN, loss.FUN = loss.FUN)
}
riskPair <- function(error1, error0, FUN = default.riskPairFUN, loss.FUN)
{
	new.riskPair(risk1 = FUN(error1), risk0 = FUN(error0), FUN = FUN, loss.FUN = loss.FUN)
}

congruity <- function(object, alternative = "greater", ...) # modified 111007
{
	if(is(object, "cvalue"))
	{
		Alt <- alt(object)
		less <- length(Alt) > 0 && Alt == "less"
		if(length(Alt) > 0 && Alt == "two.sided")
			stop("congruity incompatible with two-sided p-values")
		else if(less)
			warning(paste("alternative is", alternative))
		z <- as(object, "numeric")
		p <- pnorm(z, ...)
		pval <- if(alternative == "greater")
		{	if(less) 1 - p else p}
		else if(alternative == "two.sided")
			2 * pmin(1 - p, p)
		else
			stop("bad congruity alternative")
		names(pval) <- names(z)
		Numeric(pval)
	}
	else
		congruity(object = cvalue(object, ...), alternative = alternative)
}
confidence <- function(object, mass0 = NULL, strict.inequality = TRUE, less.than = TRUE, ...)
{
	if(is.nothing(mass0))
	{
		mass0 <- numeric(0)
	}
	else if(!is.nothing(mass0) && is(mass0, "empiricalNull"))
	{
		message("Getting local false discovery rate estimates on ", date())
		mass0 <- Mass0(mass0)#XXX|:no visible global function definition for ÔMass0Õ
	}
	else
		stopifnot(!is.nothing(mass0) && is.prob(mass0))
	con <- if(is(object, "cvalue"))
	{
		stopifnot(length(alt(object)) == 0 || alt(object) %in% c("greater", "Greater"))
		congruity(object)
	}
	stopifnot(is(mass0, "numeric") && is.prob(con))
	if(is.nothing(mass0))
	{
		if(less.than)
			con
		else
			1 - con
	}
	else if(is.prob(mass0))
	{
		assert.is(mass0, "numeric")
		if(length(mass0) == 1)
		{
			mass0 <- rep(mass0, length(con))
			names(mass0) <- names(con)
		}
		stopifnot(is.prob(mass0))
		nam <- intersect(names(mass0), names(con))
		con <- con[nam]
		mass0 <- mass0[nam]
		mass1 <- 1 - mass0
		stopifnot(length(mass1) == length(con))
		con.eq <- mass0
		con.less <- mass1 * con
		con.greater <- mass1 * (1 - con)
		stopifnot(abs(con.eq + con.less + con.greater - 1) < 1e-3)
		if(!strict.inequality)
		{
			con.less <- con.less + con.eq
			con.greater <- con.greater + con.eq
		}
		if(less.than)
			con.less
		else
			con.greater
	}
	else
		stop("cannot get observed confidence level")
}


new.alt <- function(object)
{
	if(is(object, "alt"))
		object
	else
		new("alt", object)
}
alt <- function(object)
{
	name <- "arglis"
	Alt <- "alternative"
	arglis.ok <- function(arglis){is.list(arglis) && Alt %in% names(arglis)}
	if(isS4(object) && name %in% slotNames(object))
	{
		arglis <- slot(object, name = name)
		if(!arglis.ok(arglis))
		{ message("bad arglis"); browser()}
		alt(arglis)
	}
	else if(is.list(object) && arglis.ok(object))
		alt(object[Alt][[1]])
	else if(is.character(object) && length(object) %in% c(0, 1))
		new.alt(object)
	else
		stop("cannot return alt")
}

locfdr.rname <- function(nulltype)
{
	rname <- if(nulltype == 0) # theor
		"thest"
	else if(nulltype == 1) # MLE
		"mlest"
	else if(nulltype == 2) # central matching
		"cmest"
	else
		stop("cannot compute rname for specified nulltype")
	assert.is(rname, "character")
	rname
}
setMethod("Mean", signature(x = "empiricalNull"), function(x, nulltype = default(0, "Mean nulltype"))
{
	fp0 <- s3(x)$fp0
	rname <- locfdr.rname(nulltype = nulltype)
	fp0[rname, "delta"]
})
setMethod("Sd", signature(object = "empiricalNull"), function(object, nulltype = default(0, "Sd nulltype"))
{
	fp0 <- s3(object)$fp0
	rname <- locfdr.rname(nulltype = nulltype)
	fp0[rname, "sigma"]
})

new.empiricalNull <- function(PDF, PDF0, CDF0, p0, s3, min_param = default(-Inf, "min_param"), max_param = default(Inf, "max_param"), max.p0 = 1)
{
	if(!is(CDF0, "CDF"))
		CDF0 <- new.CDF(CDF0, min_param = min_param, max_param = max_param, param.name = "z-transform of congruity", type = "probability")
	if(p0 > max.p0)
	{
		new.p0 <- max.p0
		warning(paste("p0", p0, "set to ", new.p0))
		p0 <- new.p0
	}
	if(!is(p0, "Scalar"))
		p0 <- Scalar(p0)
	new("empiricalNull", PDF = PDF, PDF0 = PDF0, CDF0 = CDF0, p0 = p0, s3 = s3)
}
silence <- function(z, nsilence, silent = NULL)
{
	assert.is(z, "numeric")
	az <- abs(z)
	get.rank <- function(x){rank(x, ties.method = "random")}
	raz <- get.rank(-az)
	maz <- max(az, na.rm = TRUE)
	if(az[raz == 1] != maz)
	{ message("bad raz"); browser()}
	stopifnot(nsilence >= 1 && floor(nsilence) == nsilence && nsilence < length(z))
	boo <- raz <= nsilence
	stopifnot(sum(boo) == nsilence)
	stopifnot(maz %in% az[boo])
	z[boo] <- if(is.null(silent))
	{
		rz <- get.rank(z)
#		oz <- order(z)
		expected.order.stats <- rnorm(ppoints(length(z)))[rz]
		stopifnot(all(length(expected.order.stats) %in% c(length(rz), length(boo))))
		new.z <- expected.order.stats[rz[boo]]
		stopifnot(length(new.z) == sum(boo))
		new.z
	}
	else if(length(silent) == 0)
		rnorm(n = sum(boo), mean = 0, sd = 1)
	else if((is.numeric(silent) || is.na(silent)) && length(silent) == 1)
		silent
	else
		stop("bad silent")
	z
}
#zaux<-try(find.package(package="locfdr", lib.loc = NULL, quiet = FALSE,verbose = FALSE),silent=T)
#if(!is.err(zaux)||length(zaux)==0){mylocfdr <- locfdr:::locfdr}
#if(is.err(zaux)||length(zaux)==0){mylocfdr <- function(zz, p.value, nulltype = 0, ...){stop("please, install package 'locfdr' available at 'cran.r-project.org/src/contrib/Archive/locfdr'");c(zz,nulltype)}}#trick 2014

mylocfdr <- function(zz, p.value, nulltype = 0, ...){stop("please, install package 'locfdr' available at 'cran.r-project.org/src/contrib/Archive/locfdr'");c(zz,nulltype)}
empiricalNull <- function(object, nulltype = default(1, "nulltype"), nsilence = 0, silent = NULL, call.browser = FALSE, cvalue.arglis = NULL, verbose = TRUE, max.p0 = 1, ...)
{
	if(is(object, "xprnSet") && is.null(cvalue.arglis))
		cvalue.arglis <- list(s3FUN = t.test)
	if(!is.null(cvalue.arglis))
	{
		assert.is(cvalue.arglis, "list")
		object <- do.call("cvalue", c(list(object, verbose = verbose), cvalue.arglis))
	}
	if(is.numeric(nulltype))
	{
		zz <- if(is(object, "cvalue"))
			as(object, "numeric")
		else if(is(object, "numeric"))
			object
		else
			stop("cannot get zz from object")
		if(nsilence > 0)
			zz <- silence(z = zz, nsilence = nsilence, silent = silent)
		boo <- is.finite(zz) # needed due to 090711 change in setAs: from = "cvalue", to = "numeric"
		if(!all(boo))
		{
			message("  empiricalNull retaining ", sum(boo), " finite z values of ", length(zz), " z values of class ", class(zz))
			zz <- zz[boo]
		}
		s3 <- try(mylocfdr(zz = zz, nulltype = nulltype, ...), silent = !verbose) # nulltype added 090711 09:19
		if(is(s3, "try-error"))
		{
			if(verbose) message("locfdr failed")
			if(call.browser)
			{
				message("press c to return NULL")
				browser()
			}
			if(verbose) message("returning NULL")
			return(NULL)
		}
		names(s3$fdr) <- names(zz)
		x <- s3$mat[, "x"]
		get.fun <- function(y)
		{
			approxfun(x, y)
		}
		PDF <- get.fun(s3$mat[, "f0theo"])
		PDF0 <- get.fun(s3$mat[, "f0"])
		rname <- locfdr.rname(nulltype = nulltype)
		fp0 <- s3$fp0
		CDF0 <- function(q)
		{
			delta <- fp0[rname, "delta"]
			sigma <- fp0[rname, "sigma"]
			pnorm(q = q, mean = delta, sd = sigma)
		}
		p0 <- fp0[rname, "p0"]
		min_param <- min(x)
		max_param <- max(x)
		stopifnot(all(is.finite(c(min_param, max_param))) && min_param <= max_param)
	}
	else
		stop("not yet implemented for non-locfdr")
	new.empiricalNull(PDF = PDF, PDF0 = PDF0, CDF0 = CDF0, p0 = p0, s3 = s3, min_param = min_param, max_param = max_param, max.p0 = max.p0)
}
assumedNull <- function(object, ...)
{
	empiricalNull(object = object, nulltype = 0, ...)
}

log_lfdr_se <- function(object, call.plot = FALSE, ...)
{
	assert.is(object, "empiricalNull")
	mat <- object@s3$mat
	fdr = mat[, "fdr"]
	lfdrse = mat[, "lfdrse"]
	fun <- approxfun(x = fdr, y = lfdrse)
	lfdr.object <- lfdr(object)
	log_lfdr_se.vec <- sapply(lfdr.object, fun)
	if(call.plot)
	{
		Mfrow()
		plot(x = fdr, y = lfdrse, ...)
		plot(lfdr.object, log_lfdr_se.vec, ...)
	}
	log_lfdr_se.vec
}
expected.lfdr <- function(object, call.plot = FALSE, ...)
{
	assert.is(object, "empiricalNull")
	lfdr.object <- lfdr(object)
	log.lfdr <- log(lfdr.object)
	sigma <- log_lfdr_se(object, call.plot = call.plot)
	elfdr <- exp(log.lfdr + sigma ^ 2 / 2)
	if(call.plot)
	{
		plot(lfdr.object, elfdr)
	}
	elfdr
}

# odds in data.s as of 100806
odds.factor <- function(object, ...)
{
	assert.is(object, "empiricalNull")
	odds(object, ...) / odds(1 - object@p0)
}

halfAlphaI <- function(x, y, ...) # see Renyi (PT, p. 595)
{
	assert.is(x, "empiricalNull")
	if(missing(y)) # relevance of empirical null
	{
		nulltype1 <- 0
		nulltype2 <- 1
		get.mean <- function(...){Mean(x, ...)}
		get.sd <- function(...){Sd(x, ...)}
		mean1 <- get.mean(nulltype = nulltype1)
		sd1 <- get.sd(nulltype = nulltype1)
		mean2 <- get.mean(nulltype = nulltype2)
		sd2 <- get.sd(nulltype = nulltype2)
	}
	else if(is(y, "empiricalNull")) # for non-ancillarity of empirical null
	{
		nulltype <- 1
		get.mean <- function(...){Mean(..., nulltype = nulltype)}
		get.sd <- function(...){Sd(..., nulltype = nulltype)}
		mean1 <- get.mean(x)
		sd1 <- get.sd(x)
		mean2 <- get.mean(y)
		sd2 <- get.sd(y)
	}
	else
		stop("bad y")
	halfAlphaRelativeEntropy(mean1 = mean1, mean2 = mean2, sd1 = sd1, sd2 = sd2)
}
new.conditioningMerit <- function(nonancillarity, relevance, nsilence)
{
	new("conditioningMerit", nonancillarity = Numeric(nonancillarity), relevance = Scalar(relevance), nsilence = Numeric(nsilence))
}
conditioningMerit <- function(object, nsilence = seq(0, length(object) / 2, 500), plot = 0, ...) # 2 ^ (0:floor(log2(length(object) / 2)))
{
	assert.is(object, "cvalue")
	assert.is(nsilence, "numeric")
	stopifnot(length(nsilence) >= 1)
	stopifnot(all(nsilence < length(object)))
	message("beginning conditioningMerit on ", date(), ".")
	get.en <- function(nsil)
	{
		empiricalNull(object, nsilence = nsil, plot = plot, ...)
	}
	en0 <- get.en(nsil = 0)
	relevance <- halfAlphaI(en0)
	nonancillarity <- sapply(nsilence, function(nsil)
	{
		en <- get.en(nsil = nsil)
		halfAlphaI(en, en0)
	})
	message("ending conditioningMerit on ", date(), ".")
	new.conditioningMerit(nonancillarity = nonancillarity, relevance = relevance, nsilence = nsilence)
}


lfdr <- function(object, zz, use.s3, factor = numeric(0), max.lfdr = Inf, ...)
{
	assert.is(object, "empiricalNull")
	if(missing(use.s3))
		use.s3 <- ("fdr" %in% names(s3(object)) && missing(zz))
	vec <- if(use.s3)
	{
		if(!missing(zz))
			warning("zz ignored by lfdr")
		s3(object)$fdr
	}
	else
	{
		warning("this conflicts with use.s3 = TRUE for an unknown reason")
		fun <- as(object, "function")
		if(missing(zz))
			fun
		else
		{
			if(is(zz, "cvalue"))
				zz <- zz@zz
			assert.is(zz, "numeric")
			sapply(zz, fun)
		}
	}
	if(length(factor) == 1)
	{
		ifelse(vec < 1, vec, exp(-abs(jitter(rep(0, length(vec)), factor = factor))))
	}
	else
		pmin(vec, max.lfdr)
}

plot_coverages <- function(seed, p0, nulltype, nsample = 200, n = 3000, sd0 = default(1, "plot_coverages sd0"))
{
	stopifnot(length(p0) == length(nulltype))
	par(mfrow = c(2, 2))
	lis <- lapply(1:length(p0), function(i)
	{
		setSeed(seed)
		cover <- coverage(nsample = nsample, n = n, p0 = p0[i], nulltype = nulltype[i], sd0 = sd0)
		plot(cover)
		cover
	})
	lis
}

new.probMass <- function(x, prob)
{
	new("probMass", x = x, prob = Numeric(prob))
}
probMass <- function(x, prob)
{
	assert.is(x, "numeric")
	stopifnot(length(x) > 0)
	assert.is(prob, "numeric")
	stopifnot(length(prob) == length(x) && all(prob > 0) && sum(prob) == 1)
	new.probMass(x = x, prob = prob)
}

new.probabilisticRisks <- function(object, file, marginal.cval)
{
	new("probabilisticRisks", object, file = file, marginal.cval = marginal.cval)
}
probabilisticRisks <- function(seed, p0, nulltype, nsample = 200, n = 3000, mean = default(2.5, "probabilisticRisks mean"), sd = default(1.25, "probabilisticRisks sd"), df = numeric(0),mean0, mean0.fixed, mean0.sd, mean0.contamination, sd0, sd0.fixed, sd0.sd, sd0.ratio, sd0.prob, df0 = numeric(0), df.prob = numeric(0), size = NULL, uncond.cval, file, verbose = FALSE, congruity.range.accept = NULL, max.abs.cvalue.error1 = Inf, sample.monitored = -1,mean0.prob,df.mass, ...) # mean and sd added shortly after 090711 15:00 to implement h090711.
{
	assert.are(list(seed, p0, nulltype, nsample, n), "numeric")
	stopifnot(length(p0) == length(nulltype))
	if(missing(mean0))
	{
		if(verbose) message("v1") # 100304
		mean0 <- if(missing(sd0) || !missing(sd0.fixed))
			default(0, "probabilisticRisks mean0")
		else
		{
			if(verbose) message("v1b")
			stop("implementation incomplete")
			if(missing(mean0.fixed))
				mean0.fixed <- default(0, "mean0.fixed")
			if(missing(mean0.sd))
				mean0.sd <- default(1, "mean0.sd")
			if(missing(mean0.prob))#XXX|: added mean0.prob in input on April2014
				mean0.contamination <- default(0.1, "mean0.contamination")
				#before April2014: mean0.contamination <- default(contamination = 0.1, "mean0.contamination") XXX|:possible error in default, unused argument (contamination = 0.1)
			assert.are(list(mean0.fixed, mean0.sd), "numeric")
			default(random.scalar(fixed = mean0.fixed, sd = mean0.sd, prob = mean0.prob), "mean0")
		}
	}
	else
		stopifnot(missing(mean0.fixed) && missing(mean0.sd))
	sd0.mass <- NULL # possibly to be reset
	if(missing(sd0))
	{
		if(verbose) message("v2") # 100304
		if(missing(sd0.fixed))
			sd0.fixed <- default(1, "sd0.fixed")
		sd0 <- if(!missing(sd0.sd) && !missing(sd0.prob) && missing(sd0.ratio))
		{
			if(verbose) message("v2a") # 100304
			stopifnot(all(sd0.sd > 0))
			stopifnot(sd0.fixed %in% sd0.sd)
			sd0.mass <- probMass(x = sd0.sd, prob = sd0.prob)
			default(random.sd(fixed = sd0.fixed, FUN = as(sd0.mass, "function")), "mixture sd0")
		}
		else if(missing(sd0.ratio))
		{
			if(verbose) message("v2b")
			if(missing(sd0.sd))
				sd0.sd <- default(1, "sd0.sd")
			assert.are(list(sd0.fixed, sd0.sd), "numeric")
			default(random.sd(fixed = sd0.fixed, sd = sd0.sd), "lognormal sd0")
		}
		else if(is.numeric(sd0.ratio) && ((length(sd0.ratio == 1) && sd0.ratio >= 1) || length(sd0.ratio) == 2) && missing(sd0.sd))
		{
			if(verbose) message("v2c")
			default(random.sd(fixed = sd0.fixed, ratio = sd0.ratio), "loguniform sd0")
		}
		else
			stop("cannot assign sd0 based on args supplied")
	}
	else
		stopifnot(missing(sd0.fixed) && missing(sd0.sd))
	if(length(df) >= 2) # added 100304: indirect computation of p-values
	{
		if(length(df.prob) == 0)
			df.prob <- rep(1 / length(df), length(df))
	}
	df.mass <- if(length(df.prob) >= 1)
	{
		stopifnot(length(df) == length(df.prob))
		probMass(x = df, prob = df.prob)
		df <- default(random.df(FUN = as(df.mass, "function")), "mixture df")#XXX|: added df.mass in input on April2014
	}
	else
		NULL
	assert.is(df0, "numeric")
	hyphen <- function(x){paste(signif(x, 2), collapse = "-")}
	sd0.text <- if(!missing(sd0.prob) && !missing(sd0.sd))
		paste("x", hyphen(sd0.sd), " prob", hyphen(sd0.prob), sep = "")
	else if(!missing(sd0.ratio))
		paste("sd0.ratio", if(length(sd0.ratio) == 1) sd0.ratio else paste(sd0.ratio, collapse = "-"), sep = "-")
	else if(!missing(sd0.ratio))
		paste(signif(sd0.fixed, 2), signif(sd0.sd, 2), sep = " ")
	else if(!missing(sd0))
		sd0
	else
		"other-sd0"
	size.unspecified <- missing(size) || is.null(size)
	if(missing(uncond.cval) || is.null(uncond.cval))
	{
		uncond.cval <- if(size.unspecified)
			NULL
		else
		{
			setSeed(seed + 1) # added 090719 22:13
			revert <- TRUE # reverts to the working method of congruity100324a.s since the method introduced in congruity100325.s fails for marginal.probabilisticRisks(size = 1e4, n = 100)
			get.cval <- function(...)
			{
				if(revert) stop("failed to revert")
				if(length(df0) == 0){df0.mass <-df.mass}
				else{stop("cannot yet assign df0.mass")}
#				df0.mass <- if(length(df0) == 0)
#					df.mass
#				else
#					stop("cannot yet assign df0.mass")
				uc.lis <- list(...)
##				if("sd0" %in% names(uc.lis))
##				{ message("redundancy of arguments?"); browser()}
				#if(verbose)
				#	message("df0.mass: ", df0.mass)
				unconditional.cvalue(size = size, mean0 = mean0, sd0.fixed = sd0.fixed, sd0.mass = sd0.mass, df0 = df0, df0.mass = df0.mass, ...) # mean = mean, sd = sd, df = df
			}#XXX|:no visible binding for global variable Ôdf0.massÕ
			if(!revert && missing(sd0.ratio))
			{
				if(verbose)
					message("trying the new method")
				get.cval(sd0 = sd0)
			}
			else
			{
				if(verbose)
					message("trying default before 24 March 2010")
				if(revert || (length(df0) == 0 && is.null(df.mass)))
				{
					message("completely reverted to default before 24 March 2010")
					unconditional.cvalue(size = size, mean0 = mean0, sd0.fixed = sd0.fixed, sd0.ratio = sd0.ratio, sd0.mass = sd0.mass)
				}
				else
				{
					#if(verbose)
					#{
					#	message("df0: ", df0.mass)
					#	message("df.mass: ", df.mass)
					#}
					get.cval(sd0.ratio = sd0.ratio) # default before 24 March 2010: unconditional.cvalue(size = size, mean0 = mean0, sd0.fixed = sd0.fixed, sd0.ratio = sd0.ratio, sd0.mass = sd0.mass)
				}
			}
		}
	}
	else if(!size.unspecified)
		stop("size and uncond.cval both specified")
	if(!is.null(uncond.cval))
	{
		main <- if(!missing(sd0.ratio))
			paste("ratio: ", sd0.ratio, sep = "")
		else if(!is.null(sd0.mass))
			paste("sd0 PMF: ", sd0.text)
		else
			sd0.mass
		plot(uncond.cval, main = main)
	}
	risk.lis <- lapply(1:length(p0), function(i)
	{
		setSeed(seed)
		probabilisticRisk(nsample = nsample, n = n, p0 = p0[i], nulltype = nulltype[i], mean = mean, sd = sd, df = df, mean0 = mean0, sd0 = sd0, df0 = df0, verbose = verbose, uncond.cval = uncond.cval, congruity.range.accept = congruity.range.accept, max.abs.cvalue.error1 = max.abs.cvalue.error1, sample.monitored = sample.monitored, ...)
	})
	file <- if(missing(file))
	{
		if(any(sample.monitored > 0))
			character(0)
		else
		{
			alt.text <- function(sca)
			{
				if(is(sca, "Scalar.factor"))
					paste("S.f", sca, sep = "-")
				else if(is.numeric(sca) && length(sca) == 1)
					sca
				else
					"other-alt"
			}
			size.text <- if(is.null(uncond.cval))
				"nominal"
			else if(is.null(size))
				"marginal"
			else
				paste("size", size, sep = "-") # added 090712
			congruity.range.accept.text <- if(is.null(congruity.range.accept))
				"risk"
			else
				paste("cra", congruity.range.accept[1], congruity.range.accept[2], sep = "-")
			ranval <- round(runif(n = 1) * 100, 0)
			fil <- try(paste("risks", seed, collapse(p0), collapse(nulltype), nsample, n, alt.text(mean), alt.text(sd), mean0, sd0.text, size.text, congruity.range.accept.text, Sys.Date(), ranval, sep = " "), silent = TRUE)
			if(is(fil, "try-error"))
				fil <- paste("risks", seed, collapse(p0), collapse(nulltype), nsample, n, "nonstandard", size.text, congruity.range.accept.text, Sys.Date(), ranval, sep = " ")
			fil
		}
	} # end file assignment
	assert.is(file, "character")
	use.file <- length(file) == 1
	risks <- new.probabilisticRisks(risk.lis, file = file, marginal.cval = uncond.cval)
	if(use.file)
	{
		saveo(risks, file = paste(file, "RData", sep = "."))
#		pdf(file = paste(file, "pdf", sep = "."))
		plot(risks, file = paste(file, "pdf", sep = "."))
	}
#	for(i in 1:length(p0))
#		plot(risks[[i]])
#	if(use.file)
#		dev.off()
	invisible(risks)
}
test.probabilisticRisks <- function(p0 = c(.95, .95, .9, .9), nulltype = c(1, 0, 1, 0), nsample = 200, n = 3000, mean = Scalar.factor(2.5), sd = Scalar.factor(1.25), ...)
{
	probabilisticRisks(seed = 34, nsample = nsample, n = n, p0 = p0, nulltype = nulltype, mean = mean, sd = sd, congruity.range.accept = c(0.01, 0.99), ...)
}
indicator.probabilisticRisks <- function(p0 = c(.95, .95), nulltype = c(1, 0), nsample = 200, ...)
{
	probabilisticRisks(seed = 34, nsample = nsample, p0 = p0, nulltype = nulltype, mean = Scalar.factor(2.5), sd = Scalar.factor(1.25), congruity.range.accept = NULL, ...)
}
stable.probabilisticRisks <- function(p0 = c(.85, .85, .9, .9, .95, .95), nulltype = c(1, 0, 1, 0, 1, 0), nsample = 500, n = stop("10000, 15000, 20000"), alpha, delta, alpha0 = alpha, delta0 = 0, nobservation, ngeneration = 1e4, ...) # ngeneration = 1e3 before 30 April 2010 # stop("ngeneration missing")
{
	test.probabilisticRisks(p0 = p0, nulltype = nulltype, nsample = nsample, n = n, alpha = alpha, delta = delta, alpha0 = alpha0, delta0 = delta0, nobservation = nobservation, ngeneration = ngeneration, ...)
}
conditional.stable.probabilisticRisks <- function(...)
{
	stable.probabilisticRisks(size = NULL, ...)
}

new.unconditional.cvalue <- function(cval, nominal.zz)
{
	new("unconditional.cvalue", cval, nominal.zz = nominal.zz)
}
unconditional.cvalue <- function(size, mean0 = default(0, "mean0"), sd0, sd0.fixed = default(1, "sd0.fixed"), sd0.ratio, sd0.mass = NULL, df0.mass = NULL, df0 = numeric(0), call.empiricalNull = FALSE,df.mass, ...) # added 090712
{
	if(!is.null(sd0.mass) && missing(sd0) && missing(sd0.ratio))
	{
		assert.is(sd0.mass, "probMass")
		sd0 <- default(random.sd(fixed = sd0.fixed, FUN = as(sd0.mass, "function")), "unconditional.cvalue mixture sd0")
	}
	else if(missing(sd0) && missing(sd0.ratio))
		stop("sd0 and sd0.ratio missing")
	else if(missing(sd0) && is.numeric(sd0.ratio) && length(sd0.ratio) %in% c(1, 2)) # sd0.ratio >= 1
	{
		sd0 <- default(random.sd(fixed = sd0.fixed, ratio = sd0.ratio), "unconditional.cvalue loguniform sd0")
	}
	else if(!(missing(sd0.ratio) && is.null(sd0.mass)))
#		stop("very better humans")
	{ message("very better humans"); browser()}
	message("\ncomputing unconditional c-values on ", date())
	if(is(sd0.mass, "probMass"))
		message("  these numeric approximations are only needed to verify exact analytic result")
	get.zz <- function(Sd)
	{
		rNorm(N = 1, Mean = value(mean0), Sd = Sd)
	}
	indices <- 1:size
	nominal.zz <- Sds <- as.numeric(rep(NA, length(indices))) # sort(sapply(indices, get.zz))
	for(i in 1:length(nominal.zz))
	{
		Sds[i] <- value(sd0)
		nominal.zz[i] <- get.zz(Sd = Sds[i])
	}
	nominal.zz <- sort(nominal.zz)
	cat("  Sds: "); print(stats(Sds))
	cat("  nominal.zz: "); print(stats(nominal.zz))
	message("finished computing unconditional c-values on ", date())
	if(any(is.na(nominal.zz)))
		stop("bad nominal.zz")
	pvalue <- indices / (size + 1)
	if(call.empiricalNull)
		en <- empiricalNull(object = nominal.zz, nulltype = 1, plot = 1, ...)
	s3FUN <- if(call.empiricalNull && missing(sd0.mass)) # s3FUN is the only important part of the return value (except for error checking)
	{
		function(zz)
		{
			qnorm(pnorm(q = zz, mean = Mean(en), sd = Sd(en)))
		}
	}
	else
	{
		true.zz <- qnorm(pvalue) # z-score that is N(0,1) under the true global null
		if(xor(is(sd0.mass, "probMass"), is(df0.mass, "probMass")))
		{
			param0.mass <- if(is(sd0.mass, "probMass") && is.null(df0.mass))
				sd0.mass
			else if(is(df.mass, "probMass") && is.null(sd0.mass))
				df0.mass##XXX|: added df.mass in input on April2014
			else
				stop("no sir")
			param0mixture <- param0.mass@x
			stopifnot(all(param0mixture > 0))
			prob.param0mixture <- param0.mass@prob
			conditional.tail.probability <- function(nom.zz) # given each component of param0mixture
			{
				stopifnot(length(nom.zz) == 1)
				get.p <- function(Sd, Df)
				{
					pNorm(Q = nom.zz, Mean = mean0, Sd = Sd, Df = Df) # pnorm(q = nom.zz, mean = mean0, sd = param0mixture) #, lower.tail = TRUE)
				}
				if(is(sd0.mass, "probMass") && is.null(df0.mass))
					get.p(Sd = param0mixture, Df = df0)
				else if(is(df.mass, "probMass") && is.null(sd0.mass))#XXX|: added df.mass in input on April2014
					get.p(Sd = sd0, Df = param0mixture)
				else
					stop("no, sir")
			}
			function(zz)
			{
				sapply(zz, function(nom.zz)
				{
					prob.higher.given.param0mixture <- conditional.tail.probability(nom.zz = nom.zz)
					if(length(prob.param0mixture) != length(prob.higher.given.param0mixture))
					{ message("bad use of param0.mass"); browser()}
					prob.higher.and.param0mixture <- prob.param0mixture * prob.higher.given.param0mixture # joint probability
					if(!all(prob.higher.and.param0mixture < prob.higher.given.param0mixture & prob.higher.and.param0mixture <= prob.param0mixture))
					{ message("bad joint probability"); browser()}
					prob.higher <- sum(prob.higher.and.param0mixture) # marginal probability
					stopifnot(prob.higher >= min(prob.higher.given.param0mixture) && prob.higher <= max(prob.higher.given.param0mixture))
					qnorm(prob.higher)
				})
			}
		}
		else
		{
			message("Without sd0.mass or df0.mass, marginal distribution will be unreliable. Press c to continue anyway.")
			browser()
			afun <- approxfun(x = nominal.zz, y = true.zz, rule = 1) # rule changed from 2 to 1 to prevent hidden errors
			function(zz)
			{
				afun(zz)
			}
		}
	} # end else for s3FUN
	stopifnot(is.function(s3FUN))
	cval <- new.cvalue(pvalue = Numeric(pvalue), s3FUN = s3FUN, arglis = list(size = size, mean0 = mean0, sd0 = sd0, alternative = "Greater"))
	new.unconditional.cvalue(cval = cval, nominal.zz = nominal.zz)
}
test.marginal.cval <- function(marginal.cval, mean, sd, sd0.fixed = 1, sd0, nulltype, seed = default(rindex(1000), "seed"), sd0.ratio, ...)
{
	assert.is(marginal.cval, "unconditional.cvalue")
	assert.are(list(nulltype), "numeric")
	if(missing(sd0) && is.numeric(sd0.ratio) && length(sd0.ratio) %in% c(1, 2)) # sd0.ratio >= 1
	{
		sd0 <- default(random.sd(fixed = sd0.fixed, ratio = sd0.ratio), "sd0")
	}
	else if(!missing(sd0.ratio))
	{ message("bad args!"); browser()}
	sub <- paste("sd0: ", if(is(sd0, "numeric")) sd0 else paste("centered at", value(sd0, random = FALSE)), "; nulltype: ", nulltype, sep = "")
	par(mfrow = c(2, 2))
#	plot(marginal.cval)
	get.cval <- function(uncond.cval)
	{
		setSeed(seed = seed)
		ra <- try(ran.cvalue(n = 3000, n0 = round(0.95 * 3000), uncond.cval = uncond.cval, mean = mean, sd = sd, sd0 = sd0, ...)) # mean = default(2.5, "mean"), sd = default(1.25, "sd"), mean0 = default(0, "mean")
		if(!is(ra, "cvalue"))
		{ message("bad ra"); browser()}
		ra
	}
	cval.uncond <- get.cval(marginal.cval)
	cval.cond <- get.cval(NULL)
	assert.are(list(cval.uncond, cval.cond), "cvalue")
	get.en <- function(cval, main)
	{
		en <- empiricalNull(object = cval, nulltype = nulltype, plot = 1)
		title(main = main, sub = sub)
		en
	}
	en.uncond <- get.en(cval.uncond, main = "marginal")
	en.cond <- get.en(cval.cond, main = "not marginal")
	assert.are(list(en.uncond, en.cond), "empiricalNull")
	get.acval <- function(cval, en)
	{
		if(!is(cval, "cvalue"))
		{ message("bad cval arg"); browser()}
		ac <- try(adjusted.cvalue(object = cval, empirical.null = en))
		if(!is(ac, "adjusted.cvalue"))
		{ message("bad ac"); browser()}
		ac
	}
	acval.uncond <- get.acval(cval = cval.uncond, en = en.uncond)
	acval.cond <- get.acval(cval = cval.cond, en = en.cond)
	assert.are(list(acval.uncond, acval.cond), "adjusted.cvalue")
	call.hist <- function(cval, xlab)
	{
		hi <- try(hist(cval@zz, xlab = xlab)) # xlim = c(-3, 3),
		if(is(hi, "try-error"))
		{ message("bad histo"); browser()}
		title(sub = sub)
	}
#	call.hist(cval = cval.cond, xlab = paste("congruity conditional on sd0 = ", if(is(sd0, "numeric")) sd0 else value(sd0, random = FALSE)))
#	call.hist(cval = cval.uncond, xlab = paste("marginal congruity"))
	plot.uncond.vs.cond <- function(cond, uncond, main)
	{
		get.vec <- function(cval){cval@zz}
		lim <- c(-3, 3)
		pl <- try(plot(get.vec(cond), get.vec(uncond), xlab = paste("congruity conditional on sd0 = ", if(is(sd0, "numeric")) sd0 else value(sd0, random = FALSE)), ylab = "marginal congruity", main = main, xlim = lim, ylim = lim, sub = sub))
		if(is(pl, "try-error"))
		{ message("bad pl"); browser()}
		abline(a = 0, b = 1, col = "gray")
	}
	plot.uncond.vs.cond(cval.cond, cval.uncond, main = paste("random"))
#	plot.uncond.vs.cond(acval.cond, acval.uncond, main = paste("adjusted by", if(nulltype == 0) "assumed" else "estimated", "null"))
	plot.adjustment <- function(random, adjusted, main)
	{
		get.vec <- function(cval){cval@zz}
		lim <- c(-3, 3)
		pl <- try(plot(get.vec(random), get.vec(adjusted), xlab = "random congruity", ylab = paste("adjusted by", if(nulltype == 0) "assumed" else "estimated", "null"), main = main, xlim = lim, ylim = lim, sub = sub))
		if(is(pl, "try-error"))
		{ message("bad pl 2"); browser()}
		abline(a = 0, b = 1, col = "gray")
	}
#	plot.adjustment(cval.uncond, acval.uncond, main = "marginal")
#	plot.adjustment(cval.cond, acval.cond, main = "not marginal")
	invisible(list(cval.uncond, cval.cond, en.uncond, en.cond, acval.uncond, acval.cond))
}


when.congruity.last.loaded <- date()

# EOF
