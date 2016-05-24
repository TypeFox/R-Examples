# binomFDR.s created by David Bickel on 17 May 2011 by moving material from Fisher.s.

BFDR <- function(alpha, P0 = 1, prob.discovery, size) # alpha = Type I test-wise error rate # ignore size
{
	stopifnot(is.prob(alpha))
	stopifnot(is.prob(P0))
	prob.false.discovery <- P0 * alpha
	stopifnot(is.prob(prob.false.discovery))
	stopifnot(is.prob(prob.discovery))
	if(all(prob.false.discovery <= prob.discovery))
		prob.false.discovery / prob.discovery
	else
		stop("contradictory arguments")
}
binom.rBFDR <- function(x, size, alpha, n, correct = default(TRUE, "binom.rBFDR correct"), ...)
{
	binom.rprob(x = x, size = size, alpha = alpha, n = n, FUN = binom.BFDR, correct = correct, ...)
}
binom.BFDR <- function(x, size, alpha, p = numeric(0), n = numeric(0), P0 = 1, max.BFDR = 1, FUN = NULL, conservative = logical(0), correct, BFDR.fun, ...)# = BFDR
# x is number of p-values <= alpha; returns p quantile of confidence posterior if p is given
{
	if(missing(correct))
		correct <- default(TRUE, "binom.BFDR correct", is.nothing(conservative))
	if("alternative" %in% names(list(...)))
		stop("alternative appeared twice in binom.BFDR")
	alternative <- if(is.nothing(conservative))
		character(0)
	else if(conservative)
		"greater"
	else
		"less"
	if(is.nothing(n) && is.null(FUN)) # deterministic BFDR estimate
	{
		prob.discovery <- if(is.nothing(p))
		{
			if(!is.nothing(conservative))
				warning(paste("conservative =", conservative, "ignored"))
			x / size
		}
		else if(is.prob(p))
			binom.prob(x = x, size = size, p = 1 - p, alternative = alternative, correct = correct, ...)
		else
			stop("bad p arg")
		prob.false.discovery <- P0 * alpha
		prob.discovery <- pmax(prob.discovery, prob.false.discovery / max.BFDR)
		BFDR.fun(alpha = alpha, P0 = P0, prob.discovery = prob.discovery, size = size) # pmin(max.BFDR, prob.false.discovery / prob.discovery)
	}
	else if(length(n) == 1 && n >= 1 && is.nothing(p)) # BFDR estimate determined by sampling the confidence posterior
	{
		if(is.nothing(FUN))
			FUN <- mean
		FUN(binom.rBFDR(x = x, size = size, alpha = alpha, n = n, conservative = conservative, correct = correct, ...))
	}
	else
		stop("bad binom.BFDR args")
}
estimated.BFDR <- function(object, alpha, nfeature, P0 = 1, p = numeric(0), n = numeric(0), ndiscovery.correction = 0, correct, verbose = FALSE, ...) # like estimated.FDR of Fisher.s
{
	if(missing(correct))
		correct <- default(TRUE, "estimated.BFDR correct", verbose = FALSE) # verbose = !("conservative" %in% names(list(...))))
	assert.is(object, "numeric")
	stopifnot(is.prob(object))
	if(missing(nfeature))
		nfeature <- default(length(object), "estimated.BFDR nfeature", verbose = length(object) < 1)
	if(missing(alpha))
		alpha <- object
	assert.is(alpha, "numeric")
	stopifnot(is.prob(alpha))
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
				message("\ncalling binom.BFDR with p == ", p, " and this ...:")
				print(list(...))
			}
			binom.BFDR(x = ndiscovery, size = nfeature, alpha = twer, p = p, n = n, P0 = P0, correct = correct, ...)
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
		new.bias.corrected.pvalue(vec, uncorrected = object)
	}
	else if(length(vec) == 1)
		Scalar(vec)
	else
	{ message("bad estimated BFDR"); browser()}
}
estimated.LFDR <- function(object, monotonic = FALSE, p = numeric(0), save.time = FALSE, verbose = FALSE, ties.method = "random", achieved.BFDR.fun = estimated.BFDR, ...)
# proven = TRUE added 110527; "proven" changed to "save.time" and default changed to FALSE on 110604
# Usage: pvals <- c(.05, .5); lfdr2 <- estimated.LFDR(); savej()
{
	assert.is(object, "numeric")
	assert.is(save.time, "logical")
	stopifnot(is.prob(object))
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
			if(is.err(rank.boo))
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
	stopifnot(is.prob(alpha))
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

default.MCP.args <- c("p", "conservative", "n", "FUN", "monotonic", "MCP.args")
default.MCP <- function(MCP = estimated.LFDR, p = numeric(0), conservative = logical(0), n = numeric(0), FUN = NULL, monotonic = logical(0), MCP.args = list())
{
	if(!is.nothing(p))
		MCP.args <- c(MCP.args, list(p = p))
	if(!is.nothing(conservative))
		MCP.args <- c(MCP.args, list(conservative = conservative))
	if(!is.nothing(n))
		MCP.args <- c(MCP.args, list(n = n))
	if(!is.nothing(FUN))
		MCP.args <- c(MCP.args, list(FUN = FUN))
	if(!is.nothing(monotonic))
		MCP.args <- c(MCP.args, list(monotonic = monotonic))
#	if(!is.nothing(MCP.args))
#	{ message("too risky! MCP.args:"); print(MCP.args); browser()} # MCP.args <- c(MCP.args, list(MCP.args = MCP.args))
	assert.is(MCP.args, "list")
	stopifnot(is.nothing(MCP.args) || (is.character(names(MCP.args)) && length(MCP.args) == length(names(MCP.args))))
	function(...)
	{
		apval <- do.call(MCP, c(list(...), MCP.args))
		assert.is(apval, "numeric")
		stopifnot(is.prob(apval))
		apval
	}
}
pvalue.posteriorP0 <- function(distr0, fam, x, individual = FALSE, unknown.param = NULL, MCP = estimated.LFDR, p = numeric(0), conservative = logical(0), n = numeric(0), FUN = NULL, monotonic = logical(0), MCP.args = list(), ...)
{
	assert.is(distr0, "univariateDistribution")
	assert.is(x, "numeric")
#	if(missing(individual) || )
#		individual <- FALSE
	get.posteriorP0 <- function(fam, individual, confidence, unknown.param)
	{
		if(is.nothing(individual) || is.nothing(confidence))
			stop("bad get.posteriorP0")
		if(is(fam, "Family") && is(x, "numeric") && length(x) >= 1 && individual)
		{ message("I really do not want to pass those arguments to posteriorP0"); browser()}
		posteriorP0(distr0 = distr0, fam = fam, x = x, individual = individual, confidence = confidence, unknown.param = unknown.param, ...)
	}
	if(is.null(unknown.param)) # estimation required
	{
		pval <- get.posteriorP0(fam = NULL, individual = TRUE, confidence = TRUE, unknown.param = NULL)
		stopifnot(length(pval) == length(x) && is.prob(pval))
		if(!is.null(MCP))
		{
			assert.is(MCP, "function")
			get.bias.corrected.pvalue <- default.MCP(MCP = MCP, p = p, conservative = conservative, monotonic = monotonic, n = n, FUN = FUN, MCP.args = MCP.args)
			assert.is(get.bias.corrected.pvalue, "function")
			apval <- as(get.bias.corrected.pvalue(pval), "numeric")
			stopifnot(length(apval) == length(pval) && is.prob(apval))
			pval <- apval
		}
		new.posteriorP0(pval, x = x, family = fam, distr0 = distr0)
	}
	else # returns true LFDR
	{
		assert.is(fam, "Family")
		get.posteriorP0(fam = fam, individual = individual, confidence = FALSE, unknown.param = unknown.param)
	}
}
pvalue.posteriorP0.Chisq <- function(df, P0, x, ...)
{
	ncp0 <- 0
	pvalue.posteriorP0(distr0 = Chisq(df = df, ncp = ncp0), fam = Family.ChisqMixture(df = df, ncp0 = ncp0), x = x, ...)
}

pvalue.posteriorP0.error <- function(distr0, fam, nsim, nfeature, unknown.param, ann = "pvalue.posteriorP0", pvalue.posteriorP0.args = list(), verbose = FALSE, ...)
{
	arglis <- list(...)
	if(verbose)
		print(arglis)
	bad.args <- default.MCP.args
	if(any(bad.args %in% names(arglis)))
	{
		message("put the following arguments inside pvalue.posteriorP0.args, not as args of pvalue.posteriorP0.error"); print(bad.args); browser()
	}
	posteriorP0.fun <- function(...)
	{
		do.call(pvalue.posteriorP0, c(list(distr0 = distr0, fam = fam, ...), pvalue.posteriorP0.args))
	}
	assert.is(posteriorP0.fun, "function")
	assert.is(unknown.param, "list")
	assert.are(list(nsim, nfeature), "numeric")
	x.fun <- function()
	{
		dis <- Distr(object = fam, unknown.param = unknown.param)
		x <- r(dis)(n = nfeature)
		assert.is(x, "r.twoDistrMixture") # required to have null.indicator slot
		x
	}
	assert.is(x.fun, "function")
	posteriorP0.error(posteriorP0.fun = posteriorP0.fun, true.posteriorP0.fun = posteriorP0.fun, nsim = nsim, nfeature = nfeature, unknown.param = unknown.param, x.fun = x.fun, null.indicator = "x", ann = ann, ...)
}
pvalue.posteriorP0.error.Chisq <- function(df, P0, ncp = ncp, nsim, nfeature, verbose = FALSE, ...)
{
	ncp0 <- 0
	arglis <- list(...)
	if(verbose)
		print(arglis)
	if(any(default.MCP.args %in% names(arglis)))
	{
		message("put the following arguments inside pvalue.posteriorP0.args, not as args of pvalue.posteriorP0.error.Chisq"); print(default.MCP.args); browser()
	}
	pvalue.posteriorP0.error(distr0 = Chisq(df = df, ncp = ncp0), fam = Family.ChisqMixture(df = df, ncp0 = ncp0), nsim = nsim, nfeature = nfeature, unknown.param = list(P0 = P0, ncp = ncp), ...)
}

pvalue.posteriorP0.errors <- function(FUN = pvalue.posteriorP0.error, ..., args.lis = lapply(pvalue.posteriorP0.args.lis, function(dummy){list()}), pvalue.posteriorP0.args.lis, name = names(pvalue.posteriorP0.args.lis), seed = 33)
{
	ok <- try(assert.are(args.lis, "list"))
	if(is.err(ok))
	{ message("bad args.lis"); browser()}
	assert.are(pvalue.posteriorP0.args.lis, "list")
	assert.is(name, "character")
	stopifnot(all(length(args.lis) == c(length(pvalue.posteriorP0.args.lis), length(name))))
	lis <- lapply(1:length(args.lis), function(i)
	{
		setSeed(seed)
		err <- do.call(FUN, c(list(..., pvalue.posteriorP0.args = pvalue.posteriorP0.args.lis[[i]]), args.lis[[i]]))
		assert.is(err, "posteriorP0.error")
		err
	})
	new.posteriorP0.errors(lis, name = name)
}
pvalue.posteriorP0.errors.Chisq <- function(...)
{
	pvalue.posteriorP0.errors(FUN = pvalue.posteriorP0.error.Chisq, ...)
}


# near end of file:

when.binomFDR.last.loaded <- date()

# EOF
