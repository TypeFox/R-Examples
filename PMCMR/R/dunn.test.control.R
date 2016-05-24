# dunn.test.control.R
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

dunn.test.control <-
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
		C <- tiesum / (12 * (n - 1))
        	return(C)
        }
        METHOD <- paste("Dunn's-test for multiple","
                         comparisons with one control", sep="\t")
        C <- getties(x.rank, n)
        if (C != 0) warning("Ties are present. z-quantiles were corrected for ties.")           
        # mean Rsum of controll is in R.bar[1]
        compare.stats <- function(j) {
            dif <- abs(R.bar[1] - R.bar[j]) 
            A <- n * (n+1) / 12
            B <- (1 / R.n[1] + 1 / R.n[j])
            zval <- dif / sqrt((A - C) * B)
            return(zval)
        }
        PSTATv <- rep(NA, k-1)
        for (j in 2:k) {PSTATv[j-1] <- compare.stats(j)}
        # unadjusted p-values
        PVALv <- 2 * pnorm(abs(PSTATv), lower.tail = FALSE)
        # adjusted p-values
        PADJv <- p.adjust(PVALv, method = p.adjust.method)

        LNAME <- levels(g)[2:k]

        # build matrix
        PSTAT <- matrix(data=PSTATv, nrow = (k-1), ncol = 1,
                        dimnames = list(LNAME, levels(g)[1]))
        PVAL <- matrix(data=PADJv, nrow = (k-1), ncol = 1,
                        dimnames = list(LNAME, levels(g)[1]))
        ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                    statistic = PSTAT, p.adjust.method = p.adjust.method)
        class(ans) <- "PMCMR"
        ans
}
