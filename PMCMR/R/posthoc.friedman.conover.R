# posthoc.friedman.conover.test.R
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

posthoc.friedman.conover.test <- function(y, ...)
    UseMethod("posthoc.friedman.conover.test")

posthoc.friedman.conover.test.default <-
    function(y, groups, blocks, p.adjust.method = p.adjust.methods, ...){
        if ((is.matrix(y)) | (is.data.frame(y))) {
            groups <- factor(c(col(y)))
            blocks <- factor(c(row(y)))
            DNAME <- paste(deparse(substitute(y)))
            GRPNAMES <- colnames(y)
        }
        else {
            if (any(is.na(groups)) || any(is.na(blocks))) 
                stop("NA's are not allowed in groups or blocks")
            if (any(diff(c(length(y), length(groups), length(blocks))))) 
                stop("y, groups and blocks must have the same length")
            if (any(table(groups, blocks) != 1)) 
                stop("Not an unreplicated complete block design")

            DNAME <- paste(deparse(substitute(y)), ",",
                           deparse(substitute(groups)), "and",
                           deparse(substitute(blocks)))
            groups <- factor(groups)
            blocks <- factor(blocks)
            GRPNAMES <- as.character(levels(groups))
        }
        p.adjust.method <- match.arg(p.adjust.method)
    	n <- length(levels(blocks))
        k <- length(levels(groups))
    	y <- y[order(groups, blocks)]
    	mat <- matrix(y, nrow = n, ncol = k, byrow = FALSE)
        for (i in 1:length(mat[, 1])) mat[i, ] <- rank(mat[i, ])
        R.sum <- colSums(mat)
        METHOD <- paste("Conover's test for a two-way","
                    balanced complete block design", sep="\t")
    
        # Friedman's T1 value
        A1 <- 0
        for (i in 1:n){
            for (j in 1:k){
                A1 <- A1 + mat[i,j]^2
            }
        }
        C1 <- (n * k * (k + 1)^2) / 4
        TT <- 0
        for (j in 1:k) {
            TT <- TT + (R.sum[j] - ((n * (k + 1))/2))^2
        }
        T1 <- ((k - 1) * TT) / (A1 - C1)
        
        A <- 2 * k * (1 - T1 / (k * (n-1))) * ( A1 - C1)
        B <- (n - 1) * (k - 1)
        # Pairwise comparisons
        compare.stats <- function(i,j) {
            dif <- abs(R.sum[i] - R.sum[j]) 
            tval <- dif / sqrt(A / B)
            return(tval)
        }
        PSTAT <- pairwise.table(compare.stats,levels(groups),
                                p.adjust.method="none" )
    
        compare.levels <- function(i,j) {
            dif <- abs(R.sum[i] - R.sum[j]) 
            tval <- dif / sqrt(A / B)
            pval <- 2 * pt(q=abs(tval), df=((n-1)*(k-1)), lower.tail=FALSE)
            return(pval)
        }
        PVAL <- pairwise.table(compare.levels,levels(groups),
                               p.adjust.method=p.adjust.method)
        
        colnames(PSTAT) <- GRPNAMES[1:(k-1)]
        rownames(PSTAT) <- GRPNAMES[2:k]
        colnames(PVAL) <- GRPNAMES[1:(k-1)]
        rownames(PVAL) <- GRPNAMES[2:k]
        ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                    statistic = PSTAT, p.adjust.method =p.adjust.method)
        class(ans) <- "PMCMR"
        ans
}
