# distribution.r created by David Bickel on 1 August 2010.

#library(distr)
##library(distrDoc) # conflicted with data.s in R version 2.9.1 (2009-06-26)
#Source("Fisher.r") # data.r") # likelihood.r")

setClass("SteinChisq", representation(df = "Scalar", ncp = "Scalar"))
setValidity("SteinChisq", function(object)
{
	ok <- length(object@df) == 1 && object@df > 0 && length(object@ncp) == 1
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})

setClass("blank.Distr", representation("NULL"))

setClass("SteinChi", representation(df = "Scalar", ncp = "Scalar"))
setValidity("SteinChi", function(object)
{
	ok <- length(object@df) == 1 && object@df > 0 && length(object@ncp) == 1
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setAs(from = "SteinChisq", to = "SteinChi", function(from)
{
	SteinChi(df = from@df, ncp = from@ncp)
})
setMethod("d", signature(object = "SteinChi"), function(object)
{
	stop("SteinChi density function not yet implemented")
})
setMethod("p", signature(object = "SteinChi"), function(object)
{
	stopifnot(object@df == 1)
	dis <- as(object, "SteinChisq")
	function(q, ...)
	{
		p(dis)(q ^ 2, ...)
	}
})
setMethod("q.r", signature(object = "SteinChi"), function(object)
{
#	stop("SteinChi quantile function not yet implemented")
	dis <- as(object, "SteinChisq")
	function(...)
	{
		sqrt(q.r(dis)(...))
	}
})
setMethod("r", signature(object = "SteinChi"), function(object)
{
	dis <- as(object, "SteinChisq")
	function(n)
	{
		ran <- r(dis)(n = n)
		if(any(is.na(ran)) || any(ran < 0))
		{ message("bad ran"); browser()}
		sqrt(ran)
	}
})
#
#setClass("SteinChisq", representation(df = "Scalar", ncp = "Scalar"))
#setValidity("SteinChisq", function(object)
#{
#	ok <- length(object@df) == 1 && object@df > 0 && length(object@ncp) == 1
#	if(!ok)
#	{ printInvalid(object); browser()}
#	ok
#})
setAs(from = "SteinChi", to = "SteinChisq", function(from)
{
	SteinChisq(df = from@df, ncp = from@ncp)
})
setMethod("d", signature(object = "SteinChisq"), function(object)
{
	stop("SteinChisq density function not yet implemented")
})
setMethod("p", signature(object = "SteinChisq"), function(object)
{
	function(q, lower.tail = TRUE, log.p = FALSE, call.browser = FALSE)
	{
		if(any(is.na(q)) || !all(q >= 0))
		{ message("bad q:"); print(stats(q)); browser()}
		legacy <- FALSE
		if(call.browser)
		{ message("SteinChisq browser"); browser()}
		if(legacy)
		{
			upperP <- pchisq(q = object@ncp, df = object@df, ncp = q, lower.tail = TRUE)
			# pcorrect <- function(n, mu.sq, x.sq){1 - pchisq(q = x.sq, df = n, ncp = mu.sq)} # from conditional-fiducial.r
			P <- if(lower.tail)
				1 - upperP
			else
				upperP
			if(log.p)
				log(P)
			else
				P
		}
		else
		{
			dis <- Chisq(df = object@df, ncp = q)
			p(dis)(q = object@ncp, lower.tail = !lower.tail, log.p = log.p)
		}
	}
})
setMethod("q.r", signature(object = "SteinChisq"), function(object)
{
	brute <- TRUE
	if(brute)
		q.SteinChisq(object = object)
	else # Venables (1975) "Calculation of Confidence Intervals for Noncentrality Parameters"
	{
		function(p, lower.tail = TRUE, log.p = FALSE)
		{
			stopifnot(is.prob(p))
			x <- object@ncp
			f <- object@df
			s <- sqrt(2 * (2 * x - f + 2)) # equation (8) of Venables (1975)
			stopifnot(all(is.finite(s)))
			get.ncp <- function(z) # equation (11) of Venables (1975)
			{
				term <- z + (z ^ 2 - 1) / s - z / s ^ 2 + ((2/3) * (f - 1) * (z ^ 2 - 1)) / s ^ 3
					+ (-(1/6) * (f - 1) * (4 * z ^ 3 - z) + (1/2) * (f - 2) * z) / s ^ 4
					+ ((4/15) * (f - 1) * (3 * z ^ 4 + 2 * z ^ 2 - 11)) / s ^ 5
					+ (-(1/90) * (f - 1) * (96 * z ^ 5 + 164 * z ^ 3 - 767 * z) - (4/9) * (f - 1) * (f - 2) * (2 * z ^ 3 - 5 * z) + (1/2) * (f - 2) * z) / s ^ 6
				x - f + 2 + s * term
			}
			z <- qnorm(p, lower.tail = lower.tail, log.p = log.p)
			stopifnot(all(is.finite(z)) && length(z) == length(p))
			qq <- get.ncp(z = z)
			stopifnot(all(is.finite(qq) & qq >= 0))
			stopifnot(length(qq) == length(p))
			qq
		}
	}
})
setMethod("r", signature(object = "SteinChisq"), function(object)
{
	r.SteinChisq(object = object)
})


setClass("absTd", representation(df = "Scalar", ncp = "scalar"))
setValidity("absTd", function(object)
{
	ok <- length(object@df) == 1 && object@df > 0 && length(object@ncp) == 1
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setAs(from = "absTd", to = "Td", function(from)
{
	Td(df = from@df, ncp = from@ncp)
})
setAs(from = "absTd", to = "Chisq", function(from)
{
	stopifnot(from@df == Inf)
	Chisq(df = 1, ncp = from@ncp ^ 2)
})
setMethod("d", signature(object = "absTd"), function(object)
{
	td <- as(object, "Td")
	function(x, log = FALSE, call.browser = FALSE)
	{
		Dfun <- function(y, take.log)
		{
			options(warn = 2)
			Den <- try(d(td)(x = y, log = take.log), silent = TRUE)
			if(is.err(Den))
			{
				Den <- if(take.log) -Inf else 0
				options(warn = 0)
				if(call.browser)
					warning(paste("bad Den set to ", Den, " for df ", object@df, " and ncp ", object@ncp, sep = ""))
			}
			options(warn = 0)
			Den
		}
		if(call.browser)
		{ message("d absTd call.browser is TRUE"); browser()}
		stopifnot(all(x >= 0))
		Dens <- Dfun(y = -x, take.log = FALSE) + Dfun(y = x, take.log = FALSE)
		zero <- Dens == 0
		log.Dens <- log(Dens)
		dens <- if(any(zero))
			ifelse(zero, pmax(Dfun(y = -x, take.log = TRUE), Dfun(y = x, take.log = TRUE)), log.Dens)
		else
			log.Dens
		if(any(dens == -Inf) && call.browser)
		{ warning("0 Likelihood!")}
		if(log)
			dens
		else
			exp(dens)
	}
})
setMethod("p", signature(object = "absTd"), function(object)
{
	if(object@ncp == 0)
	{
		td <- as(object, "Td")
		function(q, lower.tail = TRUE, log.p = FALSE, call.browser = FALSE)
		{
			if(any(is.na(q)) || !all(q >= 0))
			{ message("bad q:"); print(stats(q)); browser()}
			upperP <- 2 * pt(q = q, df = object@df, lower.tail = FALSE)
			P <- if(lower.tail)
				1 - upperP
			else
				upperP
			if(log.p)
				log(P)
			else
				P
		}
	}
	else if(object@df == Inf)
	{
		dis <- as(object, "Chisq")
		function(q, lower.tail = TRUE, log.p = FALSE, call.browser = FALSE)
		{
			if(any(is.na(q)) || !all(q >= 0))
			{ message("q failed assumptions:"); print(stats(q)); browser()}
			prob <- p(dis)(q = q ^ 2, lower.tail = lower.tail, log.p = log.p)
			if(is.na(prob))
			{ message("prob: ", prob); browser()}
			prob
		}
	}
	else
		stop("not implemented for both finite degrees of freedom and nonzero noncentrality")
})
setMethod("q.r", signature(object = "absTd"), function(object)
{
	td <- as(object, "Td")
	function(p, lower.tail = TRUE, log.p = FALSE)
	{
		stopifnot(log.p || is.prob(p))
		#if(lower.tail)
		#	stop("lower.tail = TRUE has not yet been tested")
		abs.t <- q(td)(p = p / 2, lower.tail = lower.tail, log.p = log.p)
		if(lower.tail)
			abs.t <- -abs.t
		stopifnot(all(abs.t >= 0))
		abs.t
	}
})
setMethod("r", signature(object = "absTd"), function(object)
{
	td <- as(object, "Td")
	function(n)
	{
		get.ran <- function(ran.fun){assert.is(ran.fun, "function"); abs(ran.fun(n = n))}
		ran <- if(td@param@df == Inf)
		{
			try(get.ran(ran.fun = function(n){rnorm(mean = td@param@ncp, n = n)})) # abs(rnorm(mean = td@param@ncp, n = n))
		}
		else
			get.ran(ran.fun = r(td)) # abs(r(td)(n = n))
		if(is.err(ran))
		{ message("ran Jack ran"); browser()}
		if(any(is.na(ran)) || any(ran < 0))
		{ message("bad ran"); browser()}
		ran
	}
})


setClass("absTdMixture", representation(df = "Scalar", P0 = "Scalar", ncp0 = "scalar", ncp = "scalar"))
setValidity("absTdMixture", function(object)
{
	P0.ok <- is.prob(object@P0)
	len.ok <- length(object@df) == 1 && length(object@ncp0) == 1 && length(object@ncp) == 1
	df.ok <- object@df > 0
	ncp.ok <- object@ncp0 != object@ncp
	oks <- c(P0.ok = P0.ok, len.ok = len.ok, df.ok = df.ok, ncp.ok = ncp.ok)
	ok <- all(oks)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setMethod("d", signature(object = "absTdMixture"), function(object)
{
	get.atd <- function(name){absTd(df = object@df, ncp = slot(object, name = name))}
	atd0 <- get.atd(name = "ncp0")
	atd1 <- get.atd(name = "ncp")
	function(x, log = FALSE)
	{
		stopifnot(x >= 0)
		get.Dens <- function(atd)
		{
			assert.is(atd, "absTd")
			d(atd)(x = x, log = FALSE)
		}
		P0 <- P0(object)
		P1 <- P1(object)
		Dens <- P0 * get.Dens(atd = atd0) + P1 * get.Dens(atd = atd1)
		if(log)
			log(Dens)
		else
			Dens
	}
})
setClassUnion("absTdObject", c("absTd", "absTdMixture"))
setClassUnion("univariateDistribution", c("UnivariateDistribution", "absTdObject",  "SteinChi", "SteinChisq", "blank.Distr"))#"twoDistrMixture",
setClass("twoDistrMixture", representation(P0 = "Scalar", distr0 = "univariateDistribution", distr1 = "univariateDistribution"))
setMethod("d", signature(object = "twoDistrMixture"), function(object)
{
	distr0 <- object@distr0
	distr1 <- object@distr1
	function(x, log = FALSE)
	{
		get.Dens <- function(distr)
		{
			assert.is(distr, "univariateDistribution")
			d(distr)(x = x, log = FALSE)
		}
		P0 <- P0(object)
		P1 <- P1(object)
		Dens <- P0 * get.Dens(distr = distr0) + P1 * get.Dens(distr = distr1)
		if(log)
			log(Dens)
		else
			Dens
	}
})
setMethod("p", signature(object = "twoDistrMixture"), function(object)
{
	distr0 <- object@distr0
	distr1 <- object@distr1
	max0 <- q(distr0)(p = 1)
	P0 <- object@P0
	pfun0 <- p(distr0)
	pfun1 <- p(distr1)
	distr0.le.distr1 <- pfun1(max0) == 0 # stronger than stochastic inequality
	function(q, lower.tail = TRUE, log.p = FALSE, call.browser = FALSE)
	{
		if(any(is.na(q)) || !all(q >= 0))
		{ message("bad q:"); print(stats(q)); browser()}
		stopifnot(length(q) == 1)
		lowerP <- if(q <= max0)
		{
			cond.p <- pfun0(q = q)
			cond.p * P0
		}
		else
		{
			cond.p <- 1 - pfun1(q = q)
			1 - cond.p * (1 - P0)
		}
		P <- if(lower.tail)
			lowerP
		else
			1 - lowerP
		if(log.p)
			log(P)
		else
			P
	}
})
setMethod("q.r", signature(object = "twoDistrMixture"), function(object)
{
	function(p, lower.tail = TRUE, log.p = FALSE)
	{
		if(log.p)
		{
			p <- exp(p) # stop("log.p has not yet been implemented")
			log.p <- FALSE
		}
		stopifnot(is.prob(p))
		P0 <- object@P0
		qfun0 <- q.r(object@distr0)
		qfun1 <- q.r(object@distr1)
		distr0.le.distr1 <- qfun0(p = 1) <= qfun1(p = 0) # stronger than stochastic inequality
		qq <- if(lower.tail && distr0.le.distr1)
		{
			stopifnot(length(p) == 1)
			if(p <= P0)
			{
				qfun <- qfun0
				cond.p <- p / P0
			}
			else if(1 - p <= 1 - P0)
			{
				qfun <- qfun1
				cond.p <- 1 - (1 - p) / (1 - P0)
			}
			else
				stop("a probability-zero event has occurred")
			qfun(p = cond.p, lower.tail = lower.tail, log.p = FALSE)
		}
		else
			stop("q.r, twoDistrMixture has not yet been implemented for those args")
		stopifnot(!any(is.na(qq)) && length(qq) == length(p))
		stopifnot(all(qq >= 0))
		qq
	}
})
setClass("r.twoDistrMixture", representation("numeric", null.indicator = "Numeric"))
new.r.twoDistrMixture <- function(object, null.indicator)
{
	new("r.twoDistrMixture", object, null.indicator = Numeric(null.indicator))
}
setMethod("r", signature(object = "twoDistrMixture"), function(object)
{
	function(n)
	{
		un <- runif(n = n)
		get.ran <- function(distr)
		{
			ran <- r(distr)(n = n)
			stopifnot(length(ran) == length(un))
			ran
		}
		null.boo <- un <= object@P0
		new.r.twoDistrMixture(ifelse(null.boo, get.ran(distr = object@distr0), get.ran(distr = object@distr1)), null.indicator = ifelse(null.boo, 1, 0))
	}
})

setClassUnion("abscontDistribution", c("AbscontDistribution", "absTdObject", "twoDistrMixture", "SteinChi", "SteinChisq", "blank.Distr"))
setClassUnion("discreteDistribution", c("DiscreteDistribution", "twoDistrMixture", "blank.Distr"))
setClassUnion("univariateDistribution", c("UnivariateDistribution", "absTdObject", "twoDistrMixture", "SteinChi", "SteinChisq", "blank.Distr"))

# near end of file:

#Source(file = "distribution.s") # functions

