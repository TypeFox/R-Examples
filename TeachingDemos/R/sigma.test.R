sigma.test <- 
function (x, sigma = 1, sigmasq = sigma^2,
        alternative = c("two.sided", "less", "greater"),
        conf.level = 0.95, ...) {
  alternative <- match.arg(alternative)
  sigma <- sqrt(sigmasq)
  n <- length(x)
  xs <- var(x)*(n-1)/sigma^2
  out <- list(statistic = c("X-squared" = xs))
  class(out) <- "htest"
  out$parameter <- c(df = n-1)
  minxs <- min(c(xs, 1/xs))
  maxxs <- max(c(xs, 1/xs))
  PVAL <- pchisq(xs, df = n - 1)
  
  out$p.value <- switch(alternative,
                        two.sided = 2*min(PVAL, 1 - PVAL),
                        less = PVAL,
                        greater = 1 - PVAL)
  out$conf.int <- switch(alternative,
                         two.sided = xs * sigma^2 *
1/c(qchisq(1-(1-conf.level)/2, df = n-1), qchisq((1-conf.level)/2, df
= n-1)),
                         less = c(0, xs * sigma^2 /
qchisq(1-conf.level, df = n-1)),
                         greater = c(xs * sigma^2 /
qchisq(conf.level, df = n-1), Inf))
  attr(out$conf.int, "conf.level") <- conf.level
  out$estimate <- c("var of x" = var(x))
  out$null.value <- c(variance = sigma^2)
  out$alternative <- alternative
  out$method <- "One sample Chi-squared test for variance"
  out$data.name <- deparse(substitute(x))
  names(out$estimate) <- paste("var of", out$data.name)
  return(out)}
