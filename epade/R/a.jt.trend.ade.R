a.jt.trend.ade <-
function (x, g, alternative = c("two.sided", "increasing", "decreasing"))
{
    if (!is.numeric(x))  stop("x is not numeric")
    if (!is.numeric(g) & !is.ordered(g)) stop("group should be numeric or ordered factor")
    
    alternative <- match.arg(alternative)
    n <- length(x)
    if (length(g) != n) stop("lengths of x and group don't match")
    
    ties <- length(unique(x)) != n
    #if (ties)  warning("TIES: p-value based on normal approximation")
    
# Calculation
    gsize <- table(g)
    ng <- length(gsize)
    cgsize <- c(0, cumsum(gsize))
    x <- x[order(g)]
    jtrsum <- jtmean <- jtvar <- 0
    for (i in 1:(ng - 1)) {
        na <- gsize[i]
        nb <- n - cgsize[i + 1]
        jtrsum <- jtrsum + sum(rank(x[(cgsize[i] + 1):n])[1:na]) -
            na * (na + 1)/2
        jtmean <- jtmean + na * nb/2
        jtvar <- jtvar + na * nb * (na + nb + 1)/12
    }
    jtrsum <- 2 * jtmean - jtrsum
    statistic <- jtrsum
    names(statistic) <- "JT"
    
    if (ties | (n > 50)) {
        zstat <- (statistic - jtmean)/sqrt(jtvar)
        p <- pnorm(zstat)
        p <- switch(alternative, two.sided = 2 * min(p,
            1 - p), increasing = 1 - p, decreasing = p)
    }
    else {
stop('N must be >= 50')
    }
    out <- list(statistic, alternative, as.numeric(p))
    names(out)<-c('statistic', 'alternative', 'p-Value')
    return(out)
}
