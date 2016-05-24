fisher.pval <- function (k1, n1, k2, n2,
                        alternative=c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  
  if (any(k1 < 0) || any(k1 > n1) || any(n1 <= 0)) stop("k1 and n1 must be integers with 0 <= k1 <= n1")
  if (any(k2 < 0) || any(k2 > n2) || any(n2 <= 0)) stop("k2 and n2 must be integers with 0 <= k2 <= n2")
  if (any(k1 + k2 <= 0)) stop("either k1 or k2 must be non-zero")

  l <- max(length(k1), length(n1), length(k2), length(n2)) # ensure that all vectors have the same length
  if (length(k1) < l) k1 <- rep(k1, length.out=l)
  if (length(n1) < l) n1 <- rep(n1, length.out=l)
  if (length(k2) < l) k2 <- rep(k2, length.out=l)
  if (length(n2) < l) n2 <- rep(n2, length.out=l)

  k <- k1 + k2

  pval <- switch(alternative,
                 greater   = phyper(k1 - 1, n1, n2, k, lower.tail=FALSE),
                 less      = phyper(k1, n1, n2, k, lower.tail=TRUE),
                 two.sided = 2 * pmin(phyper(k1 - 1, n1, n2, k, lower.tail=FALSE), phyper(k1, n1, n2, k, lower.tail=TRUE)))
  pval <- pmax(0, pmin(1, pval))        # clamp p-value to range [0,1] (may be > 1 in two-sided approximation)
  pval
}
