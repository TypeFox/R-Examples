cpma <- function(ps) {
    ps <- na.omit(ps)
    log.ps <- -log(ps)
    rate <- 1 / mean(log.ps)
    p.obs <- sum(dexp(log.ps, rate=rate, log=TRUE))
    p.exp <- sum(dexp(log.ps, rate=1, log=TRUE))
    chisq <- -2 * (p.obs - p.exp) * sign(rate - 1)
    structure(list(
        statistic=structure(chisq, names="X-squared"),
        parameter=structure(rate, names="lambda"),
        p.value=pchisq(chisq, df=1, lower.tail=FALSE),
        method="Cross Phenotype Meta-Analysis (CPMA)",
        data.name=deparse(substitute(ps))
    ), class="htest")
}
