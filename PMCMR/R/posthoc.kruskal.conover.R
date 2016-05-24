# posthoc.kruskal.conover.test.R
# Part of the R package: PMCMR
#
# Copyright (C) 2015, 2016 Thorsten Pohlert
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

posthoc.kruskal.conover.test <- function(x, ...) UseMethod("posthoc.kruskal.conover.test")

posthoc.kruskal.conover.test.default <-
function(x, g, p.adjust.method = p.adjust.methods, ...){
        ## taken from stats::kruskal.test
        
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
        p.adjust.method <- match.arg(p.adjust.method)
        x.rank <- rank(x)
        R.bar <- tapply(x.rank, g, mean,na.rm=T)
        R.n <- tapply(!is.na(x), g, length)
        g.unique <- unique(g)
        k <- length(g.unique)
        n <- sum(R.n)
        getties <- function(x, n) {
           x.sorted <- sort(x)
       	   pos <- 1
           tiesum <- 0
	     while (pos <= n) {
	       	val <- x.sorted[pos]  
			nt <- length(!is.na(x.sorted[x.sorted==val]))
			pos <- pos + nt
			 if (nt > 1){
			     tiesum <- tiesum + nt^3  - nt
			 }	     
                }
           	C <- 1 - tiesum / (n^3 - n)
        	C <- min(c(1,C))
        	return(C)
        }
        METHOD <- paste("Conover's-test for multiple","
                         comparisons of independent samples", sep="\t")
        C <- getties(x.rank, n)
        if (C != 1) warning("Ties are present. Quantiles were corrected for ties.")
        # Kruskal-Wallis statistic
        H <- (12 / (n * (n + 1))) * sum(tapply(x.rank, g, "sum")^2 / R.n) - 3 * (n + 1)
        H.cor <- H / C

        if (C == 1) {
            S2 <- n * (n + 1) / 12
        } else {
            S2 <-   ( 1 / (n - 1)) * (sum(x.rank^2) - (n * (((n + 1)^2) / 4)))
        }
        compare.stats <- function(i,j) {
            dif <- abs(R.bar[i] - R.bar[j]) 
            B <- (1 / R.n[i] + 1 / R.n[j])
            D <- (n - 1 - H.cor) / (n - k)
            tval <- dif / sqrt(S2 * B * D)
            return(tval)
        }
        PSTAT <- pairwise.table(compare.stats,levels(g),
                                p.adjust.method="none" )
        compare.levels <- function(i,j) {
            dif <- abs(R.bar[i] - R.bar[j]) 
            B <- (1 / R.n[i] + 1 / R.n[j])
            D <- (n - 1 - H.cor) / (n - k)
            tval <- dif / sqrt(S2 * B * D)
#            pval <- pt(abs(tval), df=n - k, lower.tail = FALSE)
            pval <- 2 * pt(abs(tval), df=n - k, lower.tail=FALSE)
            return(pval)
        }
        PVAL <- pairwise.table(compare.levels,levels(g),
                               p.adjust.method=p.adjust.method )
        
        ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                    statistic = PSTAT, p.adjust.method = p.adjust.method)
        class(ans) <- "PMCMR"
        ans
}

posthoc.kruskal.conover.test.formula <-
function(formula, data, subset, na.action,
         p.adjust.method = p.adjust.methods, ...)
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
    y <- do.call("posthoc.kruskal.conover.test", c(as.list(mf),
                                                   p.adjust.method))
    y$data.name <- DNAME
    y
}
