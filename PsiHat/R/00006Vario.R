

wilkinson.test <- function(x, mu = 0, alternative = "greater")
{
	assert.is(x, "numeric")
	p.value <- if(alternative == "greater")
	{
		prob <- function(q, lower.tail)
		{
			pr <- pt(q = q, df = length(x) - 1, lower.tail = lower.tail)
			stopifnot(length(pr) == 1)
			pr
		}
		q <- abs(mean(x)) / Sd(x, se.of.mean = TRUE)
		stopifnot(length(q) == 1 && q >= 0)
		sum(c(prob(q = q, lower.tail = FALSE), prob(q = -q, lower.tail = TRUE))) # pchisq(q = sum(x ^ 2), df = length(x), lower.tail = FALSE)
	}
	else if(alternative == "less")
		1 - wilkinson.test(x = x, mu = mu, alternative = "greater")$p.value
	else
		stop(paste("alternative '", alternative, "' not recognized", sep = ""))
	list(p.value = p.value)
}

##-------
setClass("bias.corrected.pvalue", representation("numeric", uncorrected = "Numeric", ranks = "Numeric")) # "Numeric" caused names problem
setValidity("bias.corrected.pvalue", function(object)
{
	prob.ok <- 1#is_prob(object@.Data) && is_prob(object@uncorrected)
	nam.ok <- sameNames(object, object@uncorrected)
	ranks.ok <- is_nothing(object@ranks) || (sameNames(object, object@ranks) && all(1 <= object@ranks & object@ranks <= length(object)))
	oks <- c(prob.ok = prob.ok, nam.ok = nam.ok, ranks.ok = ranks.ok)
	ok <- all(oks==T)
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})
setAs(from = "bias.corrected.pvalue", to = "numeric", function(from)
{
	vec <- from@.Data
	names(vec) <- names(from)
	vec
})
setMethod("[", signature(x = "bias.corrected.pvalue", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@.Data <- as(x, "numeric")[i]
	x@uncorrected <- x@uncorrected[i]
  x
})
new_bias.corrected.pvalue <- function(object, uncorrected, ranks = numeric(0))
{
	if(!sameNames(object, uncorrected))
	{ message("bad names 3"); browser()}
#	object <- new_Numeric(object) # "Numeric" caused names problem
	uncorrected <- nNumeric(uncorrected)
	if(!sameNames(object, uncorrected))
	{ message("bad names 4"); browser()}
#	print(names(object)); print(names(uncorrected))
	new("bias.corrected.pvalue", object, uncorrected = uncorrected, ranks = nNumeric(ranks))
}

Pbinom <- function(x, size, prob, lower.tail = TRUE, correct = default(TRUE, "Pbinom correct"), correction = if(correct) 1/2 else 0, inclusive = TRUE, verbose = FALSE)
{
	stopifnot(length(x) == 1 && length(size) == 1 && length(prob) == 1)
	stopifnot(x >= 0 && size >= 1 && x <= size) # added 090916
	get.xvec <- function(from, to){Seq(from = from, to = to, return.na.on.error = TRUE)}
	xvec <- if(lower.tail)
		get.xvec(0, if(inclusive) x else x - 1)
	else if(x <= size)
		get.xvec(if(inclusive) x else x + 1, size)
	else
		stop("bad Pbinom args")
	xvec.ok <- !any(is.na(xvec)) && all(xvec %in% 0:size)
	mass <- function(x){dbinom(x = x, size = size, prob = prob)}
	get.delta.mass <- function()
	{
		stopifnot(correct)
		correction * mass(x = x)
	}
	if(!xvec.ok && !inclusive)
	{
		Mass <- 0
		p <- if(correct)
			Mass + get.delta.mass()
		else
			Mass
		if(verbose && correct && p == 0)
		{ message("cannot correct"); browser()}
		p
	}
	else if(xvec.ok) # begin nonzero mass
	{
		Mass <- sum(sapply(xvec, mass))
		if(correct)
		{
			delta.mass <- get.delta.mass()
			corrected.p <- if(inclusive)
				Mass - delta.mass
			else
			{
				if(correction == 1/2)
					warning("arguments canceled each other's effect")
				Mass + delta.mass
			}
			if(verbose && corrected.p == 0)
			{ message("corrected.p is ", corrected.p); browser()}
			corrected.p
		}
		else
			Mass
	} # end nonzero mass
	else
		stop("bad xvec.ok")
} # end Pbinom
binom_limit <- function(x, size, p, correct = default(TRUE, "binom_limit correct"), ...) # inverse Pbinom for fixed x
{
	stopifnot(is_prob(p))
	stopifnot(length(p) == 1)
	cdf <- function(prob)
	{
		Pbinom(x = x, size = size, prob = prob, correct = correct, ...)
	}
	probs <- binom_prob(x = x, size = size, p = p, alternative = 2)
	stopifnot(length(probs) >= 2)
	stopifnot(is_prob(probs))
	opt.prob <- optimize(f = function(prob){abs(cdf(prob) - p)}, lower = min(probs), upper = max(probs))$minimum
	stopifnot(is_prob(opt.prob))
	opt.prob
}
binom_rprob <- function(x, size, n, FUN = binom_prob, ...){
	sapply(runif(n = n), function(p){FUN(x = x, size = size, p = p, ...)})
}
binom_prob <- function(x, size, p, alternative = character(0), correct = TRUE, ...)
{
	stopifnot(is_prob(p))
	if(length(p) > 1)
		return(sapply(p, function(component){binom_prob(x = x, size = size, p = component, alternative = alternative, correct = correct, ...)})) # correction = correction, call.qbinom = call.qbinom,
	stopifnot(length(p) == 1)
	if(all(p %in% c(0, 1)))
		return(p)
	if(TRUE) # missing(correction))
	{
		get.prob <- function(alternative)
		{
			conf <- if(alternative == "less")
				p
			else if(alternative == "greater")
				1 - p
			else
				stop("very bad alternative")
			i <- if(alternative == "less")
				2
			else if(alternative == "greater")
				1
			else
				stop("very bad alternative")
			limit <- sapply(conf, function(conf.level){binom.test(x = x, n = size, alternative = alternative, conf.level = conf.level, ...)$conf.int[i]})
			if(alternative == "less")
				limit
			else if(alternative == "greater")
				limit
			else
				stop("very bad alternative")
		}
		if(is_nothing(alternative))
		{
			binom_limit(x = x, size = size, p = p, correct = correct, ...)
		}
		else if(is.numeric(alternative) && alternative == 2)
			c(less = get.prob("less"), greater = get.prob("greater")) # mean
		else
			get.prob(alternative)
	}
#	else
#	if(is.numeric(correction) && length(correction) == 1 && is_nothing(alternative))
#	{
#		if(missing(correction))
#			correction <- 1 / 2
#		correction <- nscalar(correction)
#		stopifnot(is_prob(correction))
#		prob <- 1 - binom_limit(x = x, size = size, p = p, correction = correction, lower.tail = TRUE, ...) # confidence.lyx
#		if(!is_prob(prob))
#		{ message("bad prob:"); print(prob); browser()}
#		prob
#	}
	else
		stop("bad binom_prob arguments")
} # end binom_prob
monotonic.pvalue <- function(object, corrected, uncorrected, ranks = numeric(0), monotonic = TRUE)
{
	if(missing(corrected) && missing(uncorrected) && missing(ranks))
	{
		assert.is(object, "bias.corrected.pvalue")
		corrected <- as(object, "numeric")
		uncorrected <- object@uncorrected
		ranks <- object@ranks
	}
	else
		stopifnot(missing(object))
	if(is_nothing(ranks))
		ranks <- rank(uncorrected)
	if(is.null(names(corrected)))
		names(corrected) <- make.names(1:length(corrected))
	nam <- names(corrected)
	if(is.null(names(uncorrected)))
		names(uncorrected) <- nam
	if(is.null(names(ranks)))
		names(ranks) <- nam
	stopifnot(sameNames(corrected, uncorrected, ranks))
	ord <- order(ranks)
	reorder <- function(x){stopifnot(length(x) == length(ord)); x[ord]}
	sorted.uncorrected <- reorder(uncorrected)
	sorted.corrected <- reorder(corrected)
	ok <- is.character(names(sorted.corrected)) && sameNames(sorted.corrected, sorted.uncorrected)
	if(!ok)
	{ message("monotonicity problem"); browser()}
	if(is_nothing(monotonic))
	{
		for(i in length(sorted.corrected):2)
		{
			stopifnot(sorted.uncorrected[i - 1] <= sorted.uncorrected[i])
			if(sorted.corrected[i - 1] > sorted.corrected[i])
				sorted.corrected[i - 1] <- sorted.corrected[i]
		}
	}
	else if(monotonic && length(sorted.corrected) >= 2)
	{
		for(i in 2:length(sorted.corrected))
		{
			ok <- sorted.uncorrected[i - 1] <= sorted.uncorrected[i]
			if(is.na(ok) || !ok)
			{ message("sorted.uncorrected[i - 1] <= sorted.uncorrected[i] is FALSE:"); print(c(sorted.uncorrected[i - 1], sorted.uncorrected[i])); browser()}
			if(sorted.corrected[i - 1] > sorted.corrected[i])
				sorted.corrected[i] <- sorted.corrected[i - 1]
		}
	}
#	else
#		stop("that did not just happen")
	corrected <- sorted.corrected[nam]
	if(is.null(names(corrected)) || is.null(names(as(corrected, "numeric"))))
	{ message("bad names 5"); browser()}
	new_bias.corrected.pvalue(corrected, uncorrected = uncorrected, ranks = ranks)
}
#---------------------
BFDR <- function(alpha, P0 = 1, prob.discovery, size) {# alpha = Type I test-wise error rate # ignore size

	stopifnot(is_prob(alpha))
	stopifnot(is_prob(P0))
	prob.false.discovery <- P0 * alpha
	stopifnot(is_prob(prob.false.discovery))
	stopifnot(is_prob(prob.discovery))
	if(all(prob.false.discovery <= prob.discovery))
		prob.false.discovery / prob.discovery
	else
		stop("contradictory arguments")
}
binom_rBFDR <- function(x, size, alpha, n, correct = default(TRUE, "binom_rBFDR correct"), ...){

	binom_rprob(x = x, size = size, alpha = alpha, n = n, FUN = binom_BFDR, correct = correct, ...)
}
binom_BFDR <- function(x, size, alpha, p = numeric(0), n = numeric(0), P0 = 1, max.BFDR = 1, FUN = NULL, conservative = logical(0), correct, BFDR.fun, ...){# = BFDR
# x is number of p-values <= alpha; returns p quantile of confidence posterior if p is given

	if(missing(correct))
		correct <- default(TRUE, "binom_BFDR correct", is_nothing(conservative))
	if("alternative" %in% names(list(...)))
		stop("alternative appeared twice in binom_BFDR")
	alternative <- if(is_nothing(conservative))
		character(0)
	else if(conservative)
		"greater"
	else
		"less"
	if(is_nothing(n) && is.null(FUN)) # deterministic BFDR estimate
	{
		prob.discovery <- if(is_nothing(p))
		{
			if(!is_nothing(conservative))
				warning(paste("conservative =", conservative, "ignored"))
			x / size
		}
		else if(is_prob(p))
			binom_prob(x = x, size = size, p = 1 - p, alternative = alternative, correct = correct, ...)
		else
			stop("bad p arg")
		prob.false.discovery <- P0 * alpha
		prob.discovery <- pmax(prob.discovery, prob.false.discovery / max.BFDR)
		BFDR.fun(alpha = alpha, P0 = P0, prob.discovery = prob.discovery, size = size) # pmin(max.BFDR, prob.false.discovery / prob.discovery)
	}
	else if(length(n) == 1 && n >= 1 && is_nothing(p)) # BFDR estimate determined by sampling the confidence posterior
	{
		if(is_nothing(FUN))
			FUN <- mean
		FUN(binom_rBFDR(x = x, size = size, alpha = alpha, n = n, conservative = conservative, correct = correct, ...))
	}
	else
		stop("bad binom_BFDR args")
}

estimated.BFDR <- function(object, alpha, nfeature, P0 = 1, p = numeric(0), n = numeric(0),
                           ndiscovery.correction = 0, correct, verbose = FALSE, ...) {# like estimated.FDR of Fisher.s
	if(missing(correct))
		correct <- default(TRUE, "estimated.BFDR correct", verbose = FALSE) # verbose = !("conservative" %in% names(list(...))))
	assert.is(object, "numeric")
	stopifnot(is_prob(object))
	if(missing(nfeature))
		nfeature <- default(length(object), "estimated.BFDR nfeature", verbose = length(object) < 1)
	if(missing(alpha))
		alpha <- object
	assert.is(alpha, "numeric")
	stopifnot(is_prob(alpha))
	stopifnot(length(alpha) %in% c(1, length(object)))
	vec <- sapply(1:length(alpha), function(i)
	{
		twer <- alpha[i]
		ndiscovery <- pmax(0, pmin(nfeature, ndiscovery.correction + sum(object <= twer)))
		bfdr.hat <- if(twer == 0)
			1
		else
		{
			if(verbose)
			{
				message("\ncalling binom_BFDR with p == ", p, " and this ...:")
				print(list(...))
			}
			binom_BFDR(x = ndiscovery, size = nfeature, alpha = twer, p = p, n = n, P0 = P0, correct = correct, ...)
		}
		if(any(is.nan(bfdr.hat)))
		{ message("bad bfdr.hat"); browser()}
		bfdr.hat
	})
	if(any(is.nan(vec)))
	{ message("bad BFDR hat"); browser()}
	if(length(vec) == length(object))
	{
		names(vec) <- names(object)
		new_bias.corrected.pvalue(vec, uncorrected = object)
	}
	else if(length(vec) == 1)
		nscalar(vec)
	else
	{ message("bad estimated BFDR"); browser()}
}

estimated.LFDR <- function(object, monotonic = FALSE, p = numeric(0), save.time = FALSE, verbose = FALSE, ties.method = "random", achieved.BFDR.fun = estimated.BFDR, ...){
# proven = TRUE added 110527; "proven" changed to "save.time" and default changed to FALSE on 110604
# Usage: pvals <- c(.05, .5); lfdr2 <- estimated.LFDR(); savej()

	assert.is(object, "numeric")
	assert.is(save.time, "logical")
	stopifnot(is_prob(object))
#	ord <- order(object)
#	sorted.p <- object[ord]
	if(save.time)
	{
		alpha <- pmin(1, 2 * object)
	}
	else
	{
		r <- rank(object, ties.method = ties.method)
		stopifnot(!any(duplicated(r)) && all(r == floor(r)))
		alpha <- as.numeric(rep(NA, length(object)))
		for(i in 1:length(r))
		{
			R <- 2 * r[i] # - 1 corrected 110524
			rank.boo <- try(R <= length(object))
			if(is_err(rank.boo))
			{ message("bad rank.boo"); browser()}
			alpha[i] <- if(rank.boo)
			{
				boo <- R == r
				stopifnot(sum(boo) == 1 && length(boo) == length(object))
				object[boo]
			}
			else
				0
		}
	}
	stopifnot(all(is.finite(alpha)))
	stopifnot(length(alpha) == length(object))
	names(alpha) <- names(object)
	stopifnot(is_prob(alpha))
	if(verbose)
	{
		message("\ncalling achieved.BFDR.fun with save.time == ", save.time, " and this ...:")
		print(list(...))
	}
	lfdr <- achieved.BFDR.fun(object = object, alpha = alpha, verbose = verbose, p = p, ...)
	if(monotonic)
	{
		assert.is(lfdr, "bias.corrected.pvalue")
		lfdr <- monotonic.pvalue(lfdr)
	}
	lfdr
}

#---------------------
SFDR <- function(alpha, P0 , size, prob.discovery){ 

	stopifnot(is_prob(alpha))
	stopifnot(is_prob(P0))
	prob.false.discovery <- P0 * alpha
	stopifnot(is_prob(prob.false.discovery))
	stopifnot(is_prob(prob.discovery))
	if(all(prob.false.discovery <= prob.discovery))
		min(prob.false.discovery / (prob.discovery*(1-(1-alpha)^size)),1)
	else
		stop("contradictory arguments")
}



b.notrobust<-function(object,P0,...){
	stopifnot(is_prob(P0))
	estimated.BFDR (object=as(object,'Numeric'),BFDR.fun=BFDR,P0 = P0,...)}

b.robust<-function(object,P0,...){
	stopifnot(is_prob(P0))
	estimated.BFDR (object=as(object,'Numeric'),BFDR.fun=SFDR,P0 = P0,...)}

#--------------------

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

setClass("CDF", representation("function", min_param = "scalar", max_param = "scalar", param.name =  "character", type = "character")) # Cf. "pposterior" of logical.r
setValidity("CDF", function(object)
{
	ok <- object@min_param <= object@max_param && length(object@param.name) == 1 && length(object@type) >= 1
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
	alternative <- nalt(object)
	alt.ok <- is(alternative, "alt") && (length(alternative) == 0 || alternative != "two.sided")
	cla.ok <- !is(object@zz, "Numeric")
	ok <- alt.ok && cla.ok
	if(!ok)
	{ printInvalid(object); browser()}
	ok
})

#
setMethod("[", signature(x = "cvalue", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x@zz <- x@zz[i]
	stopifnot(validObject(x))
	x
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



#
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
	new_cvalue(pvalue = from@pvalue, s3FUN = t.test, arglis = list(alternative = from@alternative))
})




#
new_ttest <- function(pvalue, stat, df, alternative, level1, level2)
{
	new("ttest", pvalue = nNumeric(pvalue), stat = stat, df = nNumeric(df), alternative = alternative, level1 = level1, level2 = level2)
}
new_empiricalNull <- function(PDF, PDF0, CDF0, p0, s3, min_param = default(-Inf, "min_param"), max_param = default(Inf, "max_param"), max.p0 = 1)
{
	if(!is(CDF0, "CDF"))
		CDF0 <- new_CDF(CDF0, min_param = min_param, max_param = max_param, param.name = "z-transform of congruity", type = "probability")
	if(p0 > max.p0)
	{
		new_p0 <- max.p0
		warning(paste("p0", p0, "set to ", new_p0))
		p0 <- new_p0
	}
	if(!is(p0, "Scalar"))
		p0 <- nScalar(p0)
	new("empiricalNull", PDF = PDF, PDF0 = PDF0, CDF0 = CDF0, p0 = p0, s3 = s3)
}
new_cvalue <- function(pvalue, zz, s3FUN, arglis)
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
ncvalue <- function(object, s3FUN, alternative, ttest.arglis, verbose = TRUE, ...) # s3FUN, e.g., t.test, has arguments x, alternative, and possibly y; returns list with "pvalue" in names
{
	if(missing(ttest.arglis))
	{
		if(verbose)
			message("\ncomputing cvalue object on ", date())
		if(missing(alternative))
			alternative <- default("Greater", "\n  cvalue alternative", verbose = verbose)
		alternative <- nalt(alternative)
		if(missing(s3FUN) || !is.function(s3FUN))
			stop(paste('s3FUN missing or wrong class; supply argument s3FUN, e.g., t.test,', 'which takes arguments x, alternative, and possibly y; returns list with "pvalue" in names'))
		pvalue <- PValue(object = object, FUN = s3FUN, alternative = alternative, verbose = verbose, ...)
		if(!is.na(pvalue[1]) && all(pvalue == pvalue[1], na.rm = TRUE))
		{ message("all p-values the same"); browser()}
		new_cvalue(pvalue = pvalue, s3FUN = s3FUN, arglis = list(alternative = alternative, ...))
	}
	else if(missing(s3FUN) && missing(alternative))
	{
		assert.is(ttest.arglis, "list")
		t.obj <- do.call("nttest", c(list(object), ttest.arglis))
		as(t.obj, "cvalue")
	}
	else
		stop("cvalue arg err")
}
mylocfdr <- locfdr#function(zz, p.value, nulltype = 0, ...){stop("please, install package 'locfdr' available at 'cran.r-project.org/src/contrib/Archive/locfdr'");c(zz,nulltype)}
nempiricalNull <- function(object, nulltype = default(1, "nulltype"), nsilence = 0, silent = NULL, call.browser = FALSE, cvalue.arglis = NULL, verbose = TRUE, max.p0 = 1, ...)
{
	if(is(object, "nxprnSet") && is.null(cvalue.arglis))
		cvalue.arglis <- list(s3FUN = t.test)
	if(!is.null(cvalue.arglis))
	{
		assert.is(cvalue.arglis, "list")
		object <- do.call("ncvalue", c(list(object, verbose = verbose), cvalue.arglis))
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
	new_empiricalNull(PDF = PDF, PDF0 = PDF0, CDF0 = CDF0, p0 = p0, s3 = s3, min_param = min_param, max_param = max_param, max.p0 = max.p0)
}
assumedNull <- function(object, ...)
{
	nempiricalNull(object = object, nulltype = 0, ...)
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


new_CDF <- function(object, min_param, max_param, param.name, type)
{
	new("CDF", object, min_param = nscalar(min_param), max_param = nscalar(max_param), param.name = param.name, type = type)
}
nCDF <- function(type, ...)
{
	fun <- if(type == "confidence")
		confidence_CDF
	else if(type == "probability")
		probability_CDF
	else
		stop("bad type of CDF")
	fun(...)
}

blank_CDF <- function(object, param.name = "no parameter")
{
	if(missing(object))
		object <- "no function specified for CDF"
	assert.is(object, "character")
	fun <- function(param)
	{
		stop(object)
	}
	new_CDF(object = fun, min_param = -Inf, max_param = Inf, param.name = param.name, type = "no type")
}


PValue <- function(object, get.PValue, alternative = default("Greater", "alternative"), verbose = TRUE, ...)
{
	if(missing(get.PValue))
		get.PValue <- PValueFUN(alternative = alternative, ...)
	recall <- function(object){PValue(object = object, get.PValue = get.PValue, alternative = alternative, ...)}
	if(is(object, "numeric"))
	{
		pval <- get.PValue(object)
		names(pval) <- names(object)
	}
	else if(is(object, "matrix"))
	{
		if(verbose)
			message("\n  began p-value computations ", date())
		pval <- apply(object, 1, get.PValue)
		if(verbose)
			message("\  completed p-value computations ", date(), "\n")
		stopifnot(nrow(object) == length(pval))
		names(pval) <- rownames(object)
	}
	#else if(is(object, "nxprnSetPair"))
	#{
	#	get.mat <- function(es)
	#	{
	#		if(is(es, "NxprnSet"))
	#			es <- logb(es)
	#		nexprs(es)
	#	}
	#	xmat <- nexprs(object@x)
	#	ymat <- nexprs(object@y)
	#	stopifnot(all(rownames(xmat) == rownames(ymat)))
	#	nam <- rownames(xmat)
	#	if(verbose)
	#		message("\n  began p-value computations (2-sample test) ", date())
	#	pval <- sapply(1:length(nam), function(i)
	#	{
	#		get.vec <- function(mat)
	#		{
	#			vec <- try(as.numeric(mat[i, ]))
	#			if(is_err(vec))
	#			{ message("vec problem"); browser()}
	#			vec
	#		}
	#		get.PValue(x = get.vec(xmat), y = get.vec(ymat))
	#	})
	#	if(verbose)
	#		message("\  completed p-value computations ", date(), "\n")
	#	stopifnot(length(nam) == length(pval))
	#	names(pval) <- nam
	#}
	#else if(is(object, "NxprnSet"))
	#	pval <- recall(object = logb(object))
	#else if(is(object, "nxprnSet"))
	#{
	#	pval <- recall(object = nexprs(object))
	#	nam <- nfeatureNames(object)
	#	stopifnot(length(pval) == length(nam))
	#	names(pval) <- nam
	#}
	else
		stop("bad class for PValue")
	if(!is.na(pval[1]) && all(pval == pval[1], na.rm = TRUE))
	{ message("all p-values same"); browser()}
	nNumeric(pval)
}

PValueFUN <- function(FUN, alternative, ...)
{
	arglis <- list(...)
	null.arg.name <- "mu"
	mu <- if(null.arg.name %in% names(arglis))
		arglis[null.arg.name][[1]]
	else
		0
	stopifnot(is.numeric(mu) && length(mu) == 1)
  if(missing(FUN))
    FUN <- t.test
  if(!is.function(FUN))
  { message("bad FUN"); browser()}
  if(missing(alternative))
  	alternative <- character(0)
  function(x, y)
  {
  	if(missing(y))
  		y <- NULL
  	insignificant.p <- if(length(alternative) == 0)
     	0.5 # 1 before 20 August 2007
    else
     	1
  	get.p <- function(alternative)
  	{
  		extract.p <- function(...){FUN(alternative = alternative, ...)$p.value}
  		if(length(alternative) == 0 || alternative == "Greater") # one-sided, but conservative toward 0.5 rather than toward 1
  		{
  			less <- get.p(alternative = "less")
  			greater <- get.p(alternative = "greater")
  			if(length(alternative) == 0)
  			{
					if(less < greater)
						less
					else if(greater < less)
						1 - greater
					else
						insignificant.p
				}
				else if(alternative == "Greater") # added 090716; corrects some conservatism
					mean(c(1 - less, greater))
				else
					stop("bad alternative")
  		}
  		else if(is.null(y))
  			extract.p(x = x)
  		else
  			extract.p(x = x, y = y)
  	}
  	effective.sample.size <- function()
  	{
  		ss <- sample_size # function(vec){sum(!is.na(vec))}; generalized 12 May 2008
  		if(is.null(y))
  			ss(x)
  		else if(identical(FUN, wilcox.test) && min(ss(x), ss(y)) >= 1)
  			max(ss(x), ss(y))
  		else
  			min(ss(x), ss(y))
  	}
    p <- if(effective.sample.size() >= 2)
		{
			present.part <- function(vec){vec[!is.na(vec)]}
			x <- present.part(x)
			if(!is.null(y)) y <- present.part(y)
			is_essentially.constant <- function(vec){all(vec[1] == vec)}
			t.test.like <- identical(FUN, t.test) || identical(FUN, wilkinson.test)
			if(t.test.like && (is_essentially.constant(x) && (is.null(y) || is_essentially.constant(y))))
			{
				sample.mean <- if(is.null(y))
					x[1]
				else
					x[1] - y[1]
				if(sample.mean == mu)
					insignificant.p
				else
				{
					significant.p <- if(sample.mean > mu && length(alternative) == 0) 1 else 0 # 0 before 20 August 2007
					warning(paste("returning p-value of ", significant.p, " for ", paste(x, collapse = " | ")), " since t.test fails.", sep = "")
					significant.p
				}
			}
			else
				get.p(alternative = alternative) # typical computation of the p-value
		}
		else
			as.numeric(NA) # insignificant.p before 070828 15:57
    nscalar(p)
  } # end function(x, y)

new_ttest <- function(pvalue, stat, df, alternative, level1, level2)
{
	new("ttest", pvalue = nNumeric(pvalue), stat = stat, df = nNumeric(df), alternative = alternative, level1 = level1, level2 = level2)
}
nttest <- function(x, y, factor.name, level1, level2, alternative = "greater", ...)
{
	if(is(x, "nXprnSet"))
		x <- logb(x)
	if(!missing(y) && is(y, "nXprnSet"))
		y <- logb(y)
	assert.is(x, "nxprnSet")
	if(missing(y) && !missing(factor.name))
	{
		get.sub <- function(level)
		{
			nxprnSubset(object = x, level = level, factor.name = factor.name)
		}
		stopifnot(length(level2) == 1)
		get.ttest <- function(lev1)
		{
			stopifnot(length(lev1) == 1)
			nttest(x = get.sub(level = lev1), y = get.sub(level = level2), level1 = lev1, level2 = level2, alternative = alternative, ...)
		}
		if(length(level1) == 1)
			get.ttest(lev1 = level1)
		else if(length(level1) == 2)
		{
			t1 <- get.ttest(lev1 = level1[1])
			t2 <- get.ttest(lev1 = level1[2])
			com <- Combine(t1, t2, xlab = level1[1], ylab = level1[2])
			if(length(com) != length(t1) + length(t2))
			{ message("length(com) != length(t1) + length(t2)"); browser()}
			com
		}
		else
			stop("bad length(level1)")
	}
	else
	{
		nam <- nfeatureNames(x)
		x <- nexprs(x)
		y <- if(missing(y))
		{
			NULL
		}
		else if(missing(factor.name))
		{
			assert.is(y, "nxprnSet")
			stopifnot(all(nam == nfeatureNames(y)))
			nexprs(y)
		}
		else
			stop("bad args!!!!!!!!!!!!!")
		if(is.null(y))
		{
			level1 <- "x"
			level2 <- "missing"
		}
		else
			stopifnot(nrow(x) == nrow(y))
		get.s3 <- function(i)
		{
			test <- function(...)
			{
				t.test(x = x[i, ], alternative = alternative, ...)
			}
			if(is.null(y))
				test(...)
			else
				test(..., y = y[i, ])
		}
		pvalue <- stat <- df <- as.numeric(rep(NA, nrow(x)))
		for(i in 1:length(stat))
		{
			s3 <- get.s3(i)
			pvalue[i] <- s3$p.value
			stat[i] <- s3$statistic
			df[i] <- s3$parameter
		}
		stopifnot(length(pvalue) == length(nam))
		names(pvalue) <- names(stat) <- names(df) <- nam
		new_ttest(pvalue = pvalue, stat = stat, df = df, alternative = alternative, level1 = level1, level2 = level2)
	}
}

#-------

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
		get_pvalue <- function(...)
		{
			t.test(..., mu = param, alternative = "greater")$p.value
		}
		if(nsample == 2)
			get_pvalue(x = object[[1]], y = object[[2]])
		else if(nsample == 1)
			get_pvalue(x = object)
		else
			stop(paste("cannot compute", param.name, "p-value from given object of class", class(object)))
	}
	confidence_CDF(object = object, pvalue.fun = pvalue.fun, param.name = param.name, min_param = -Inf, max_param = Inf, ...)
}
confidence_CDF <- function(object, pvalue.fun, param.name, min_param, max_param, ...) # cf. CD.cdf.binom from logical.r; s3.test.fun removed
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
		{ message("confidence_CDF missing object"); browser()}
		cdf <- function(param) # integrate needs vector param
		{
#			if(length(param) != 1)
#			{ message("bad param length in confidence_CDF"); browser()}
			confid <- try(ifelse(param < min_param, 0, ifelse(param > max_param, 1, pvalue.fun(object = object, param = param, ...)))) # as of 090509, no longer <= (for integrate)
			if(is(confid, "try-error"))
			{ message("Does pvalue.fun allow vector param?"); print(list(object = object, param = param, ...)); browser()}
			confid
		}
		new_CDF(object = cdf, min_param = min_param, max_param = max_param, param.name = param.name, type = "confidence")
	}
}
probability_CDF <- function(object, param.name = default("parameter", "param.name"), min_param, max_param, ...)
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
	new_CDF(object = cdf, min_param = min_param, max_param = max_param, param.name = param.name, type = "probability")
}
new_alt <- function(object)
{
	if(is(object, "alt"))
		object
	else
		new("alt", object)
}
nalt <- function(object)
{
	name <- "arglis"
	Alt <- "alternative"
	arglis.ok <- function(arglis){is.list(arglis) && Alt %in% names(arglis)}
	if(isS4(object) && name %in% slotNames(object))
	{
		arglis <- slot(object, name = name)
		if(!arglis.ok(arglis))
		{ message("bad arglis"); browser()}
		nalt(arglis)
	}
	else if(is.list(object) && arglis.ok(object))
		nalt(object[Alt][[1]])
	else if(is.character(object) && length(object) %in% c(0, 1))
		new_alt(object)
	else
		stop("cannot return alt")
}

CombinedNames <- function(object, xlab = "x", ylab = "y")
{
	if(is(object, "character"))
	{
		if(xlab == ylab)
			stop("CombinedNames xlab == ylab")
		dup0 <- duplicated(object)
		if(any(dup0))
		{ message("duplication in object:"); print(dup0); browser()}
		subnam <- function(lab)
		{
			stopifnot(length(lab) == 1)
			paste(lab, object, sep = ".")
		}
		nam <- c(subnam(xlab), subnam(ylab))
		stopifnot(length(nam) == 2 * length(object))
		dup <- duplicated(nam)
		if(any(dup))
		{ message("duplication in nam:"); print(dup); browser()}
		nam
	}
	else if(is.null(object))
	{
		warning("No names were found; returning the object")
		NULL
	}
	else # adds combined names to object
		stop("bad class")
#	{
#		names(object) <- CombinedNames(object = names(object), xlab = xlab, ylab = ylab)
#		object
#	}
}

#removeMethods("Combine")Combine <- function(x, y, ...){printGeneric("Combine"); browser()}
setGeneric("Combine", function(x,y,...) standardGeneric("Combine"))
setMethod("Combine", signature(x = "Vector", y = "Vector"), function(x, y, ...)
{
	combo <- c(x, y, ...)
	if(sameNames(x, y))
		names(combo) <- CombinedNames(names(x))
	if(is(names(combo), "character"))
		names(combo) <- make.names(names(combo), unique = TRUE)
	combo
})
setMethod("Combine", signature(x = "Numeric", y = "Numeric"), function(x, y, ...)
{
	num <- function(object){as(object, "numeric")}
	nNumeric(Combine(x = num(x), y = num(y), ...))
})
setMethod("Combine", signature(x = "ttest", y = "ttest"), function(x, y, xlab = stop("no xlab"), ylab = stop("no ylab"), ...) # like x = "cvalue", y = "cvalue"
{
	assert.is(xlab, "character")
	assert.is(ylab, "character")
	stopifnot(length(xlab) == 1)
	stopifnot(length(ylab) == 1)
	if(xlab == ylab)
		stop("Combine xlab == ylab")
	same <- function(name)
	{
		identical(slot(x, name), slot(y, name))
	}
	stopifnot(same("alternative"))
	stopifnot(same("level2"))
	len <- length(x) + length(y)
	nam <- if(all(names(x) == names(y)))
		CombinedNames(names(x), xlab = xlab, ylab = ylab)
	else
		stop("x and y have different names")
	get.vec <- function(name)
	{
		get.subvec <- function(object, lab)
		{
			z <- try(slot(object, name = name))
			if(is(z, "try-error"))
			{ message("slot(object, name = name) error!"); browser()}
			names(z) <- paste(lab, names(z), sep = ".") # sapply(names(z), function(na){paste(lab, na, sep = ".")})
			z
		}
		p1 <- get.subvec(x, lab = xlab)
		p2 <- get.subvec(y, lab = ylab)
		zvalue <- c(p1, p2)
		dup <- duplicated(names(zvalue))
		if(any(dup))
		{ message("duplication in names(zvalue):"); print(dup); browser()}
		if(any(nam != names(zvalue)))
		{ message("failed CombinedNames"); browser()}
		zvalue
	}
	x@pvalue <- nNumeric(get.vec(name = "pvalue"))
	x@stat <- get.vec(name = "stat")
	x@df <- nNumeric(get.vec(name = "df"))
	x@level1 <- paste(x@level1, y@level1, sep = ".")
	stopifnot(validObject(x))
	if(length(x) != len)
	{ message("length(x) != len"); browser()}
	x
})



#removeMethods("Rank")Rank <- function(object, ...){printGeneric("Rank"); browser()}
setGeneric("Rank", function(object,...) standardGeneric("Rank"))
setMethod("Rank", signature(object = "numeric"), function(object, factor, ...)
{
	if(missing(factor))
	{
		factor <- "warn"
		message("Rank factor = ", factor)
	}
	if(factor == "warn")
	{
		recurse <- function(factor){Rank(object = object, factor = factor, ...)}
		ra <- recurse(factor = FALSE) # no jitter
		if(all(ra == floor(ra)))
			ra
		else
		{
			warning(paste("ties broken randomly on ", date()))
			recurse(factor = 1)
		}
	}
	else
	{
		if(!is.logical(factor) || factor)
			object <- jitter(object, factor = factor)
		ra <- rank(object, ...)
		stopifnot(all(ra >= 0, na.rm = TRUE))
		nNumeric(ra)
	}
})
#setMethod("Rank", signature(object = "biasEstimate"), function(object)
#{
#	if(!object@sorted)
#		stop(paste(class(object), "object is not sorted"))
#	if(length(object@Rank) == 0)
#	{
#		vec <- nNumeric(1:length(object))
#		names(vec) <- names(object)
#		nNumeric(vec)
#	}
#	else
#		object@Rank
#})


#------------------------------------------------------------------
sorted <- function(object, ...)
{
	Slot(object = object, name = "sorted", ...)
}
#removeMethods("sorted")
#removeMethods("sample_size")sample_size <- function(object, ...){printGeneric("sample_size"); browser()}

setGeneric("sample_size", function(object,...) standardGeneric("sample_size"))

setMethod("sample_size", signature(object = "numeric"), function(object)
{
	sum(!is.na(object))
})

setMethod("sample_size", signature(object = "matrix"), function(object)
{
	size <- sapply(1:nrow(object), function(i){sample_size(as.numeric(object[i, ]))})
	names(size) <- rownames(object)
	size
})
setMethod("sample_size", signature(object = "nxprnSet"), function(object)
{
	mat <- nexprs(object)
	size <- sapply(1:nrow(mat), function(i){sample_size(as.numeric(mat[i, ]))})
	names(size) <- nfeatureNames(object)
	size
})

#--------
}