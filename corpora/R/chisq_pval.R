chisq.pval <- function (k1, n1, k2, n2, correct=TRUE,
                        alternative=c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  
  if (alternative == "two.sided") {
    X2 <- chisq(k1, n1, k2, n2, correct=correct)
    pval <- pchisq(X2, df=1, lower.tail=FALSE)
  } else {
    z <- chisq(k1, n1, k2, n2, correct=correct, one.sided=TRUE)
    pval <- pnorm(z, lower.tail = (alternative == "less"))
  }

  pval
}
