z.score.pval <- function (k, n, p=.5, correct=TRUE, alternative=c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  
  z <- z.score(k, n, p=p, correct=correct)
  pval <- switch(alternative,
                 two.sided = 2 * pnorm(abs(z), lower.tail=FALSE),
                 less      = pnorm(z, lower.tail=TRUE),
                 greater   = pnorm(z, lower.tail=FALSE))
  pval
}
