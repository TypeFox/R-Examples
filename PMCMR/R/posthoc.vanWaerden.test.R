#  posthoc.vanWaerden.test.R
#
#  Copyright (C) 2015, 2016 Thorsten Pohlert
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

posthoc.vanWaerden.test <- function(x, ...) UseMethod("posthoc.vanWaerden.test")

posthoc.vanWaerden.test.default <-
function(x, g,  p.adjust.method = p.adjust.methods, ...)
{
    if (is.list(x)) {
        if (length(x) < 2L)
            stop("'x' must be a list with at least 2 elements")
        DNAME <- deparse(substitute(x))
        x <- lapply(x, function(u) u <- u[complete.cases(u)])
        k <- length(x)
        l <- sapply(x, "length")
        if (any(l == 0))
            stop("all groups must contain data")
        g <- factor(rep(1 : k, l))
        x <- unlist(x)
    }
    else {
        if (length(x) != length(g))
            stop("'x' and 'g' must have the same length")
        DNAME <- paste(deparse(substitute(x)), "and",
                       deparse(substitute(g)))
        OK <- complete.cases(x, g)
        x <- x[OK]
        g <- g[OK]
        if (!all(is.finite(g)))
            stop("all group levels must be finite")
        g <- factor(g)
        k <- nlevels(g)
        if (k < 2)
            stop("all observations are in the same group")
    }

    n <- length(x)
    if (n < 2)
        stop("not enough observations")
    r <- rank(x)
    p.adjust.method <- match.arg(p.adjust.method)
    # transform to z-scores
    zscores <- qnorm(r / (n+1))
    AJ <- tapply(zscores, g, sum)
    NJ <- tapply(zscores, g, length)
    s2 <- (1 / (n - 1)) * sum(zscores^2)
    STATISTIC <- (1 / s2) * sum(AJ^2 / NJ)
    PARAMETER <- k - 1
    A.mn <- AJ / NJ

    compare.stats <- function(i,j) {
        dif <- abs(A.mn[i] - A.mn[j]) 
        B <- (1 / NJ[i] + 1 / NJ[j])
        tval <- dif / sqrt(s2 * (n-1-STATISTIC)/(n-k) * B)
            return(tval)
        }
    PSTAT <- pairwise.table(compare.stats,levels(g),
                                p.adjust.method="none" )
    
    compare.levels <- function(i,j) {
        dif <- abs(A.mn[i] - A.mn[j]) 
        B <- (1 / NJ[i] + 1 / NJ[j])
        tval <- dif / sqrt(s2 * (n-1-STATISTIC)/(n-k) * B)
        pval <- 2 * pt(abs(tval), df=n - k, lower.tail=FALSE)
        return(pval)
    }
    PVAL <- pairwise.table(compare.levels,levels(g),
                           p.adjust.method=p.adjust.method )
    METHOD <- paste("van der Waerden normal scores test for", "
               multiple comparisons of independent samples", sep="\t")
    
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method)
    class(ans) <- "PMCMR"
    return(ans)
}

posthoc.vanWaerden.test.formula <-
function(formula, data, subset, na.action, p.adjust.method = p.adjust.methods, ...)
{
    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
                 
   if(missing(formula) || (length(formula) != 3L))
        stop("'formula' missing or incorrect")
    mf <- eval(mf, parent.frame())  
    if(length(mf) > 2L)
       stop("'formula' should be of the form response ~ group")
    DNAME <- paste(names(mf), collapse = " by ")
    p.adjust.method <- match.arg(p.adjust.method)
    names(mf) <- NULL
    y <- do.call("posthoc.vanWaerden.test", c(as.list(mf), p.adjust.method))
    y$data.name <- DNAME
    y
}
