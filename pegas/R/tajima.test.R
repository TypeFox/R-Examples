## tajima.test.R (2014-06-26)

##   Test of the Neutral Mutation Hypothesis

## Copyright 2009-2014 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

tajima.test <- function(x)
{
    n <- if (is.list(x)) length(x) else dim(x)[1]
    if (n < 4) {
        warning("Tajima test requires at least 4 sequences")
        return(list(D = NaN, Pval.normal = NaN, Pval.beta = NaN))
    }
    khat <- mean(dist.dna(x, "N"))
    S <- length(seg.sites(x))
    if (!S) {
        warning("no segregating sites")
        return(list(D = NaN, Pval.normal = NaN, Pval.beta = NaN))
    }
    tmp <- 1:(n - 1)
    a1 <- sum(1/tmp)
    a2 <- sum(1/tmp^2)
    b1 <- (n + 1)/(3 * (n - 1))
    b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
    c1 <- b1 - 1/a1
    c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)
    D <- (khat - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
    Dmin <- (2/n - 1/a1)/sqrt(e2)
    Dmax <- ((n + 1)/(2 * n) - 1/a1)/sqrt(e2)
    tmp1 <- 1 + Dmin * Dmax
    tmp2 <- Dmax - Dmin
    a <- -tmp1 * Dmax/tmp2
    b <- tmp1 * Dmin/tmp2
    p <- pbeta((D - Dmin)/tmp2, b, a)
    p <- if (p < 0.5) 2 * p else 2 * (1 - p)
    list(D = D, Pval.normal = 2 * pnorm(-abs(D)), Pval.beta = p)
}
