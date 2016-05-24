# distribution.s created by David Bickel on 1 August 2010.

new.blank.Distr <- blank.Distr <- function(){new("blank.Distr")}

SteinChi <- new.SteinChi <- function(df = 1, ncp)
{
	new("SteinChi", df = Scalar(df), ncp = Scalar(ncp))
}
SteinChisq <- new.SteinChisq <- function(df = 1, ncp)
{
	new("SteinChisq", df = Scalar(df), ncp = Scalar(ncp))
}
q.SteinChisq <- function(object, lower = 0, upper, length.out = 1e4, ...)
{
	assert.is(object, "SteinChisq")
	stopifnot(object@df == 1)
	if(missing(upper))
		upper <- max(100, length.out * object@ncp / 100) # object@df
	stopifnot(length(lower) == 1 && length(upper) == 1)
	domain0 <- seq(lower, upper, length.out = length.out)
	pfun0 <- p(object)
	pfun.lower <- function(q)
	{
		get.p <- function(qq){pfun0(q = qq, lower.tail = TRUE)}
#		min.p <- get.p(qq = 0)
		get.p(qq = q)
	}
	pfun.upper <- function(q){stop("not ready for lower.tail = FALSE"); pfun0(q = q, lower.tail = FALSE)}
	function(p, lower.tail = TRUE, log.p = FALSE)
	{
		if(log.p) stop("log.p = TRUE is not yet supported")
		pfun <- if(lower.tail)
			pfun.lower
		else
			pfun.upper
		qfun <- inversefun(object = pfun, domain0 = domain0, ...)
		if(lower.tail)
		{
			min.p <- pfun(lower)
			q0 <- qfun(p)
			stopifnot(length(p) == length(q0))
			qq <- ifelse(p <= min.p, lower, q0) # if(p <= min.p) lower else qfun(p)
			stopifnot(length(p) == length(qq))
			ok <- !any(is.na(qq)) && all(qq >= min(domain0) & qq < max(domain0))
			if(!ok)
			{ message("bad qq"); browser()}
			qq
		}
		else
			stop("lower.tail = FALSE has not yet been implemented")
	}
}
r.SteinChisq <- function(object, brute, ...)
{
	if(missing(brute))
		brute <- TRUE
	qfun <- if(brute)
		q.SteinChisq(object = object, ...)
	else if(length(list(...)) == 0)
	{
		stop("Venables (1975) is disabled")
		q.r(object)
	}
	else
		stop("bad r.SteinChisq args")
	assert.is(qfun, "function")
	function(n)
	{
		p <- runif(n = n)
		stopifnot(is.prob(p))
		ran <- qfun(p) # sapply(p, qfun)
		if(any(is.na(ran)) || any(ran < 0))
		{ message("bad ran"); browser()}
		ran
	}
}

absTd <- new.absTd <- function(df, ncp = 0)
{
	new("absTd", df = Scalar(df), ncp = scalar(ncp))
}
d.absTd <- function(x, df, ncp = 0, ...)
{
	at <- absTd(df = df, ncp = ncp)
	d(at)(x = x, ...)
}
q.absTd <- function(p, df, ncp = 0, ...) # q.absTd(p = c(.05, 0.07418), df = 2)
{
	stopifnot(is.prob(p))
	at <- absTd(df = df, ncp = ncp)
	q.r(at)(p = p, ...)
}
two.sided.p <- function(p, alternative)
{
	if(alternative == "two.sided")
		p
	else
		ifelse(p < 0.5, 2 * p, 2 * (1 - p))
}
stat.absTd <- function(p, alternative = "two.sided", lower.tail = FALSE, ...) # stat.absTd(p = c(.05, 0.07418), df = 2) # for converting p-values to absolute values of Student t statistics
{
	p2 <- two.sided.p(p = p, alternative = alternative)
	q.absTd(p2, lower.tail = lower.tail, ...)
}

new.absTdMixture <- function(df, P0, ncp0 = 0, ncp) # deprecated 110224
{
	new("absTdMixture", df = Scalar(df), P0 = Scalar(P0), ncp0 = scalar(ncp0), ncp = scalar(ncp))
}
absTdMixture <- function(df, P0, ncp0 = 0, ncp, return.twoDistrMixture = TRUE, Distr.fun = absTd)
{
	if(P0 > 1)
	{
		warning("P0 > 1 reset to 1")
		P0 <- 1
	}
	if(P0 < 0)
	{
		warning("P0 < 0 reset to 0")
		P0 <- 0
	}
	stopifnot(length(df) == 1 && length(ncp0) == 1 && length(ncp) == 1 && length(P0) == 1 && is.prob(P0))
	if(return.twoDistrMixture)
	{
		get.distr <- function(NCP)
		{
			distr <- try(Distr.fun(df = df, ncp = NCP))
			if(is.err(distr))
			{ message("bad distr; check class(NCP):", class(NCP)); browser()}
			distr
		}
		new.twoDistrMixture(P0 = P0, distr0 = get.distr(NCP = ncp0), distr1 = get.distr(NCP = ncp))
	}
	else if(identical(Distr.fun, absTd))
		new.absTdMixture(df = df, P0 = P0, ncp0 = ncp0, ncp = ncp) # default before 110224
	else
		stop("incompatible arguments; consider return.twoDistrMixture = TRUE")
}
ChisqMixture <- function(...)
{
	absTdMixture(return.twoDistrMixture = TRUE, Distr.fun = Chisq, ...)
}

new.twoDistrMixture <- function(P0, distr0, distr1)
{
	new("twoDistrMixture", P0 = Scalar(P0), distr0 = distr0, distr1 = distr1)
}
twoDistrMixture <- function(object, ...)
{
	if(is(object, "numeric") && length(object) == 1 && is.prob(object))
		new.twoDistrMixture (P0 = P0, ...)
	else
		stop("bad args for twoDistrMixture")
}

P0 <- function(object, ...)
{
	stopifnot(length(list(...)) == 0)
	object@P0
}
P1 <- function(object, ...)
{
	1 - P0(object, ...)
}
setMethod("P0", signature(object = "SteinChi"), function(object, ...)
{
	p(object)(0)
})
setMethod("P0", signature(object = "SteinChisq"), function(object, ...)
{
	p(object)(0)
})
setMethod("P0", signature(object = "twoDistrMixture"), function(object, ...)
{
	object@P0
})
setMethod("P0", signature(object = "univariateDistribution"), function(object, ...)
{
	if("P0" %in% slotNames(object))
		object@P0
	else
		0
})

has.Distr.fun <- function(object, Distr.fun)
{
	assert.is(object, "Family")
	assert.is(Distr.fun, "function")
	boo <- sameFunction(object@Distr.fun, Distr.fun)
	if(length(boo) == 0 || is.na(boo) || !is.logical(boo))
	{ message("has.Distr.fun error"); browser()}
	boo
}
is.Family.Td <- function(object)
{
	boo <- has.Distr.fun(object = object, Distr.fun = Td)
	if(!is.logical(boo))
	{ message("is.Family.Td error"); browser()}
	boo
}
is.Family.absTd <- function(object)
{
	boo <- has.Distr.fun(object = object, Distr.fun = absTd)
	if(!is.logical(boo))
	{ message("is.Family.absTd error"); browser()}
	boo
}
is.Family.absTdMixture <- function(object)
{
	boo <- has.Distr.fun(object = object, Distr.fun = absTdMixture)
	if(!is.logical(boo))
	{ message("is.Family.absTdMixture error"); browser()}
	boo
}
is.Family.absTdObject <- function(object)
{
	boo <- try(is.Family.absTd(object) || is.Family.absTdMixture(object))
	if(is.err(boo) || !is.logical(boo))
	{ message("is.Family.absTdObject error"); browser()}
	boo
}







# near end of file:

when.distribution.last.loaded <- date()

# EOF
