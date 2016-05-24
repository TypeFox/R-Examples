binom.pval <- function (k, n, p=.5, alternative=c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  if (any(k < 0) || any(k > n) || any(n < 1)) stop("arguments must be integer vectors with 0 <= k <= n")
  if (any(p < 0) || any(p > 1)) stop("null hypothesis proportion p must be in range [0,1]")
  
  pval <- switch(alternative,
                 two.sided = 2 * pmin(pbinom(k, n, p, lower.tail=TRUE), pbinom(k-1, n, p, lower.tail=FALSE)),
                 less      = pbinom(k, n, p, lower.tail=TRUE),
                 greater   = pbinom(k-1, n, p, lower.tail=FALSE))
  pval <- pmax(0, pmin(1, pval))        # clamp p-value to range [0,1] (may be > 1 in two-sided approximation)
  pval
}
