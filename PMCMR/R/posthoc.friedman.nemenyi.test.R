# posthoc.friedman.nemenyi.test.R
# Part of the R package: PMCMR
#
# Copyright (C) 2014, 2015, 2016 Thorsten Pohlert
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


posthoc.friedman.nemenyi.test <- function(y, ...) UseMethod("posthoc.friedman.nemenyi.test")

posthoc.friedman.nemenyi.test.default <-
function(y, groups, blocks, ...){
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

	 DNAME <- paste(deparse(substitute(y)), ",", deparse(substitute(groups)), "and", deparse(substitute(blocks)) )
         groups <- factor(groups)
         blocks <- factor(blocks)
         GRPNAMES <- as.character(levels(groups))
        }
    	n <- length(levels(blocks))
        k <- length(levels(groups))
    	y <- y[order(groups, blocks)]
    	mat <- matrix(y, nrow = n, ncol = k, byrow = FALSE)
    	 for (i in 1:length(mat[, 1])) mat[i, ] <- rank(mat[i, ])
        R.mnsum <- colMeans(mat)
        p.adjust.method = "none"
        METHOD <- paste("Nemenyi multiple comparison test","
             with q approximation for unreplicated blocked data", sep="\t")
         compare.stats <- function(i,j) {
            dif <- abs(R.mnsum[i] - R.mnsum[j])
            qval <- dif / sqrt(k * (k + 1) / (6 * n))
            return(qval)
        }
        PSTAT <- pairwise.table(compare.stats,levels(groups), p.adjust.method="none" ) * sqrt(2)
         ## Set df to Inf
        ###PVAL <- 1 - ptukey(PSTAT, nmeans=k, df=1000000)
        PVAL <- 1 - ptukey(PSTAT, nmeans=k, df=Inf)
        colnames(PSTAT) <- GRPNAMES[1:(k-1)]
        rownames(PSTAT) <- GRPNAMES[2:k]
        colnames(PVAL) <- GRPNAMES[1:(k-1)]
        rownames(PVAL) <- GRPNAMES[2:k]
        ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
               statistic = PSTAT, p.adjust.method = p.adjust.method)
        class(ans) <- "PMCMR"
        ans
}

posthoc.friedman.nemenyi.test.formula <-
function(formula, data, subset, na.action, ...)
{
    if(missing(formula))
        stop("formula missing")
    ## <FIXME>
    ## Maybe put this into an internal rewriteTwoWayFormula() when
    ## adding support for strata()
    if((length(formula) != 3L)
       || (length(formula[[3L]]) != 3L)
       || (formula[[3L]][[1L]] != as.name("|"))
       || (length(formula[[3L]][[2L]]) != 1L)
       || (length(formula[[3L]][[3L]]) != 1L))
        stop("incorrect specification for 'formula'")
    formula[[3L]][[1L]] <- as.name("+")
    ## </FIXME>
    m <- match.call(expand.dots = FALSE)
    m$formula <- formula
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1L]] <- quote(stats::model.frame)
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " and ")
    names(mf) <- NULL
    y <- do.call("posthoc.friedman.nemenyi.test", as.list(mf))
    y$data.name <- DNAME
    y
}
