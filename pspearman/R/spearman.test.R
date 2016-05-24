# This file was created using file src/library/stats/R/cor.test.R
# from the R package, http://www.R-project.org

spearman.freq <- function(n)
{
	m <- n*(n^2 - 1)/6 +1
	freq <- double(m)
	g <- spearman.list[[n]]
	l <- length(g)
	freq[1:l] <- g
	freq[m + 1 - (1:l)] <- g
	freq
}

pspearman <- function(s, n, lower.tail = TRUE,
		approximation = c("exact", "AS89", "t-distribution"))
{
	approximation <- match.arg(approximation)
	if (approximation == "exact") {
		freq <- spearman.freq(n)
		if (!lower.tail) {
			s <- n*(n^2 - 1)/3 - s
		}
		k <- s/2 + 1
		sum(freq[1:k])/sum(freq)
	} else if (approximation == "AS89") # call of AS 89
		.C("prho",
		   as.integer(n),
		   as.double(s + 2*lower.tail),
		   p = double(1),
		   integer(1),
		   as.logical(lower.tail),
		   PACKAGE = "pspearman")$p
	else { # approximation using t_{n-2}
		r <- 1 - 6 * s / (n*(n^2-1)) # careful for overflow
		pt(r / sqrt((1 - r^2)/(n-2)), df = n-2, lower.tail= !lower.tail)
	}
}

spearman.test <- function(x, y, alternative = c("two.sided", "less", "greater"),
		approximation = c("exact", "AS89", "t-distribution"))
{
	alternative <- match.arg(alternative)
	approximation <- match.arg(approximation)
	DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

	if (length(x) != length(y))
		stop("'x' and 'y' must have the same length")
	OK <- complete.cases(x, y)
	x <- x[OK]
	y <- y[OK]
	n <- length(x)

	PVAL <- NULL
	NVAL <- 0

	if (n < 2)
		stop("not enough finite observations")
	PARAMETER <- NULL
	TIES <- (min(length(unique(x)), length(unique(y))) < n)
	method <- "Spearman's rank correlation rho"
	names(NVAL) <- "rho"
	r <- cor(rank(x), rank(y))
	ESTIMATE <- c(rho = r)
	if (!is.finite(ESTIMATE)) {  # all x or all y the same
		ESTIMATE[] <- NA
		STATISTIC <- c(S = NA)
		PVAL <- NA
	}
	else {
		## Use the test statistic S = sum(rank(x) - rank(y))^2
		## In the case of no ties, S = (1-rho) * (n^3-n)/6.
		s <- (n^3 - n) * (1 - r) / 6
		if (!TIES) {
			s <- 2*round(s/2)
		}
		STATISTIC <- c(S = s)
		if (TIES && approximation == "exact") {
			approximation <- "AS89"
			warning("Cannot compute exact p-values with ties")
		}
		if (approximation == "exact" && n > length(spearman.list)) {
			approximation <- "AS89"
			warning("Exact values not available for n > ", length(spearman.list))
		}
		PVAL <- switch(alternative,
				"two.sided" = {
					p <- if (s > (n^3 - n) / 6)
						pspearman(s, n, lower.tail = FALSE, approximation)
						else
						pspearman(s, n, lower.tail = TRUE, approximation)
					min(2 * p, 1)
				},
				"greater" = pspearman(s, n, lower.tail = TRUE, approximation),
				"less" = pspearman(s, n, lower.tail = FALSE, approximation))
	}

	RVAL <- list(statistic = STATISTIC,
				 parameter = PARAMETER,
				 p.value = as.numeric(PVAL),
				 estimate = ESTIMATE,
				 null.value = NVAL,
				 alternative = alternative,
				 method = method,
				 data.name = DNAME)
	class(RVAL) <- "htest"
	RVAL
}

