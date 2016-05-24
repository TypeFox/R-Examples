# Fisher.s created by David Bickel on 23 April 2011.

new.bias.corrected.pvalue <- function(object, uncorrected, ranks = numeric(0))
{
	if(!sameNames(object, uncorrected))
	{ message("bad names 3"); browser()}
#	object <- Numeric(object) # "Numeric" caused names problem
	uncorrected <- Numeric(uncorrected)
	if(!sameNames(object, uncorrected))
	{ message("bad names 4"); browser()}
#	print(names(object)); print(names(uncorrected))
	new("bias.corrected.pvalue", object, uncorrected = uncorrected, ranks = Numeric(ranks))
}
bias.corrected.pvalue.cdf <- function(ranks, nfeature, strict.equality = default(FALSE, "bias.corrected.pvalue.cdf strict.equality"), ...)
{
	assert.are(list(ranks, nfeature), "numeric")
	stopifnot(length(nfeature) == 1)
	stopifnot(length(ranks) %in% c(1, nfeature))
	stopifnot(all(ranks >= 1 & ranks <= nfeature))
	function(q)
	{
		stopifnot(length(ranks) %in% c(1, length(q))) # length(q) %in% c(1, nfeature))
		stopifnot(is.prob(q))
		p.strict.equality <- function(qu, ra)
		{
			pbeta(q = qu, shape1 = ra, shape2 = nfeature - ra + 1, ...) # bias.corrected.pvalue(object = q, ...)
		}
		p.inequality <- function(qu, ra)
		{
			stopifnot(length(qu) == 1)
			stopifnot(length(ra) == 1)
			mean(sapply(1:ra, function(ra){p.strict.equality(qu = qu, ra = ra)})) # h110425a
		}
		p <- if(strict.equality)
			p.strict.equality(qu = q, ra = ranks)
		else if(length(ranks) >= 1)
		{
			len <- max(length(ranks), length(q))
			if(length(ranks) == 1)
				ranks <- rep(ranks, len)
			if(length(q) == 1)
				q <- rep(q, len)
			stopifnot(all(sapply(list(q, ranks), length) == len))
			sapply(1:len, function(k){p.inequality(qu = q[k], ra = ranks[k])})
		}
		else
			stop("bad, very very bad")
		stopifnot(is.prob(p))
		p
	}
}
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
	if(is.nothing(ranks))
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
	if(is.nothing(monotonic))
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
	new.bias.corrected.pvalue(corrected, uncorrected = uncorrected, ranks = ranks)
}
bias.corrected.pvalue <- function(object, nfeature, strict.equality = default(FALSE, "bias.corrected.pvalue strict.equality"), monotonic = TRUE, ties.method = "random", ...)
{
	assert.is(object, "numeric")
	stopifnot(is.prob(object))
	ranks <- rank(object, ties.method = ties.method, ...)
	if(missing(nfeature))
		nfeature <- default(length(object), "bias.corrected.pvalue nfeature", verbose = length(object) <= 1)
	if(is.null(names(object)))
		names(object) <- make.names(1:length(object))
	nam <- names(object)
	stopifnot(all(sapply(list(object, ranks), length) %in% c(1, nfeature)))
	vec <- bias.corrected.pvalue.cdf(ranks, nfeature = nfeature, strict.equality = strict.equality)(q = object) # pbeta(q = object, shape1 = ranks, shape2 = nfeature - ranks + 1)
	stopifnot(length(vec) == nfeature)
#	sapply(1:nfeature, function(k)
#	{
#		i <- ranks[k]
#		pbeta(q = object[k], shape1 = i, shape2 = nfeature - i + 1)
#	})
#	object@.Data <- vec
	names(vec) <- nam # names(object) <-
	if(is.nothing(monotonic) || monotonic)
	{
		vec <- as(monotonic.pvalue(corrected = vec, uncorrected = object, ranks = ranks, monotonic = monotonic), "numeric")
	} # end if(monotonic)
	stopifnot(all(names(vec) == nam))
	if(!all(names(object) == nam))
	{ message("bad names 1"); browser()}
	if(!sameNames(vec, object))
	{ message("bad names 2"); browser()}
	names(ranks) <- nam
	new.bias.corrected.pvalue(vec, uncorrected = object, ranks = ranks)
}

# BFDR functions moved to binomFDR.s on 17 May 2011

new.bias.corrected.zvalue <- function(object, uncorrected)
{
	bcz <- new("bias.corrected.zvalue", object, uncorrected = uncorrected)
	assert.is(bcz, "bias.corrected.zvalue")
	if(is(bcz, "bias.corrected.pvalue"))
		stop("no way")
	bcz
}
setMethod("zeta", signature(object = "bias.corrected.pvalue"), function(object, ...)
{
	new.bias.corrected.zvalue(zeta(as(object, "numeric"), ...), uncorrected = zeta(object@uncorrected, ...))
})

new.bias.corrected.certainty <- function(lower, upper, nonmonotonic)
{
	new("bias.corrected.certainty", lower = Numeric(lower), upper = Numeric(upper), nonmonotonic = nonmonotonic)
}
bias.corrected.certainty <- function(object, ...)
{
	if(is(object, "bias.corrected.pvalue") && is.nothing(list(...)))
	{
		get.bound <- function(monotonic)
		{
			as(monotonic.pvalue(object, monotonic = monotonic), "numeric")
		}
		new.bias.corrected.certainty(lower = get.bound(monotonic = logical(0)), upper = get.bound(monotonic = TRUE), nonmonotonic = object)
	}
	else if(is(object, "numeric") && !is(object, "bias.corrected.pvalue"))
	{
		bcp <- bias.corrected.pvalue(object = object, monotonic = FALSE, ...)
		bias.corrected.certainty(bcp)
	}
	else
		stop("bad bias.corrected.certainty args")
}

independent.Sidak <- function(object, nfeature)
{
	assert.is(object, "numeric")
	stopifnot(is.prob(object))
	if(missing(nfeature))
		nfeature <- default(length(object), "independent.Sidak nfeature", verbose = length(object) <= 1)
	vec <- 1 - (1 - object) ^ nfeature # Efron (2010, p. 36)
	names(vec) <- names(object)
	new.bias.corrected.pvalue(vec, uncorrected = object)
}
probability.nonzero.discovery <- function(ndiscovery, nfeature, exact = TRUE)
{
	prob.no.discovery <- if(exact)
	{
		stopifnot(length(nfeature) == 1)
		prob.nondiscovery <- 1 - ndiscovery / nfeature
		stopifnot(is.prob(prob.nondiscovery))
		prob.nondiscovery ^ nfeature
	}
	else
		dpois(x = 0, lambda = ndiscovery)
	stopifnot(length(prob.no.discovery) == length(ndiscovery))
	names(prob.no.discovery) <- names(ndiscovery)
	1 - prob.no.discovery
}
estimated.FDR <- function(object, alpha, nfeature, estimate.BFDR = FALSE, exact = logical(0), control = TRUE, monotonic, P0 = 1, max.monitored.ndiscovery = -Inf) # like estimated.BFDR of binomFDR.s
{
	assert.is(object, "numeric")
	stopifnot(is.prob(object))
	if(missing(nfeature))
		nfeature <- default(length(object), "estimated.FDR nfeature", verbose = length(object) <= 1)
	if(missing(alpha))
		alpha <- object
	assert.is(alpha, "numeric")
	stopifnot(is.prob(alpha))
	stopifnot(length(alpha) %in% c(1, length(object)))
	raw.vec <- sapply(1:length(alpha), function(i)
	{
		ndiscovery <- sum(object <= alpha[i])
		expected.FDP <- P0 * alpha[i] * nfeature / ndiscovery
		if(estimate.BFDR)
		{
			if(is.nothing(exact))
				exact <- TRUE
			prob.nonzero <- probability.nonzero.discovery(ndiscovery = ndiscovery, nfeature = nfeature, exact = exact)
			BFDR <- expected.FDP / prob.nonzero
			if(ndiscovery <= max.monitored.ndiscovery)
				print(round(c(ndiscovery = ndiscovery, nfeature = nfeature, prob.nonzero = prob.nonzero, expected.FDP = expected.FDP, BFDR = BFDR), 3))
			BFDR
		}
		else if(is.nothing(exact))
			expected.FDP
		else
			stop("bad args for estimated.FDR")
	})
	vec <- pmin(1, raw.vec)
	if(length(vec) == length(object))
	{
		names(vec) <- names(object)
		fdr.hat <- new.bias.corrected.pvalue(vec, uncorrected = object)
		if(control && missing(monotonic))
			fdr.hat
		else if(!control)
		{
			if(missing(monotonic))
				monotonic <- TRUE
			fdr.hat@.Data <- as(monotonic.pvalue(fdr.hat, monotonic = monotonic), "numeric")
			fdr.hat
		}
		else
			stop("bad estimated.FDR args")
	}
	else if(length(vec) == 1 && control && missing(monotonic))
		Scalar(vec)
	else
	{ message("bad estimated FDR"); browser()}
}
estimated.truth <- function(object, FDR, control = TRUE, monotonic, return.logical = FALSE, ...)
{
	stopifnot(is.prob(FDR) && length(FDR) == 1)
	FDR <- Scalar(FDR)
	fdr.hat <- estimated.FDR(object = object, alpha = object, control = control, monotonic = monotonic, ...)
	get.discovery <- function(fdr.estimate)
	{
		stopifnot(all(fdr.estimate >= 0) && length(fdr.estimate) == length(object))
		fdr.estimate <= FDR
	}
	discovery <- get.discovery(fdr.hat)
	vec <- if(return.logical)
		discovery
	else
		ifelse(discovery, FDR, 1)
	stopifnot(is.prob(vec))
	stopifnot(length(object) == length(vec))
	names(vec) <- names(object)
	vec
}



# near end of file:

when.Fisher.last.loaded <- date()

# EOF
