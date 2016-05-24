lrp_binom <-
    function(x, n, pc0=.5, pg=.5, alternative="greater")
{
    phat <- delimit(x/n, lower=pg, upper=1)
    logLikMax <- dbinom(x=x, size=n, prob=phat, log=TRUE)
    logLikNull <- dbinom(x=x, size=n, prob=pc0, log=TRUE)
    Stat <- sign(phat - pc0) * sqrt(2 * (logLikMax - logLikNull))
    normPval(Stat, alternative=alternative)
}

Waldp_binom <-
    function(x, n, pc0=.5, pg=.5, alternative="greater")
{
    phat <- delimit(x/n, lower=pg, upper=1)
    Stat <- (phat - pc0) / sqrt(phat*(1 - phat)/n)
    normPval(Stat, alternative=alternative)
}

scorep_binom <-
    function(x, n, pc0=.5, pg=.5, alternative="greater")
{
    phat <- delimit(x/n, lower=pg, upper=1)
    Stat <- (phat - pc0) / sqrt(pc0*(1 - pc0)/n)
    normPval(Stat, alternative=alternative)
}

binomPwr <-
    function(pcA, pc0=.5, n, pg=.5, alpha=.05,
             alternative=c("two.sided", "greater", "less"),
             statistic=c("likelihood", "score", "Wald"))
{
    alt <- match.arg(alternative)
    stat <- match.arg(statistic)

    tol <- 1e-10
    x <- 0:n
    d <- dbinom(x=x, size=n, prob=pcA)
    x <- x[d > tol]
    d <- d[d > tol]
    fun <- switch(stat,
                  "likelihood" = lrp_binom,
                  "score" = scorep_binom,
                  "Wald" = Waldp_binom)
    pvals <- fun(x=x, n=n, pc0=pc0, pg=pg, alternative=alt)
    ## Power:
    sum(d[pvals < alpha])
}

d.primePwr2 <-
    function(d.primeA, d.prime0 = 0, sample.size, alpha = 0.05,
             method = c("duotrio", "tetrad", "threeAFC", "twoAFC",
             "triangle"),
             alternative = c("difference", "similarity", "two.sided",
             "less", "greater"),
             statistic = c("likelihood", "score", "Wald"))
{
    ## match and test arguments:
    alt <- match.arg(alternative)
    stat <- match.arg(statistic)
    method <- match.arg(method)
    stopifnot(length(d.primeA) == 1 && is.numeric(d.primeA) &&
              d.primeA >= 0)
    stopifnot(length(d.prime0) == 1 && is.numeric(d.prime0) &&
              d.prime0 >= 0)
    stopifnot(is.numeric(sample.size) &&
              length(sample.size) == 1 &&
              sample.size > 0)
    if(max(abs(sample.size - round(sample.size))) > 1e-6)
        warning("non-integer 'sample.size' rounded to ", round(sample.size))
    size <- round(sample.size)
    stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
              alpha > 0 && alpha < 1)
    ## Test that d.prime0, d.primeA conforms with 'alternative'
    if(alt %in% c("difference", "greater") && d.primeA < d.prime0)
        stop("Need d.primeA >= d.prime0 when alternative is '", alt, "'")
    if(alt %in% c("similarity", "less") && d.primeA > d.prime0)
        stop("Need d.primeA <= d.prime0 when alternative is '", alt, "'")
    if(alt == "difference") alt <- "greater"
    if(alt == "similarity") alt <- "less"

    pA <- psyfun(d.primeA, method=method)
    p0 <- psyfun(d.prime0, method=method)
    pg <- if(method %in% c("duotrio", "twoAFC")) 1/2 else 1/3

    binomPwr(pcA=pA, pc0=p0, n=size, pg=pg, alpha=alpha,
             alternative=alt, statistic=stat)
}

normPval <-
    function(statistic, alternative = c("two.sided", "less", "greater"))
{
    alternative <- match.arg(alternative)
    stopifnot(all(!is.na(statistic)))
    p.value <-
        switch(alternative,
               "greater" = pnorm(statistic, lower.tail = FALSE),
               "less" = pnorm(statistic, lower.tail = TRUE),
               "two.sided" = 2 * pnorm(abs(statistic), lower.tail = FALSE))
    return(p.value)
}
