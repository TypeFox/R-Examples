#  posthoc.quade.test.R
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

posthoc.quade.test <- function(y, ...) UseMethod("posthoc.quade.test")

posthoc.quade.test.default <-
function(y, groups, blocks, dist=c("TDist", "Normal"),
         p.adjust.method = p.adjust.methods,  ...)
{
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
    k <- nlevels(groups)
    b <- nlevels(blocks)
    y <- matrix(unlist(split(y, blocks)), ncol = k, byrow = TRUE)
    y <- y[complete.cases(y), ]
    dist <- match.arg(dist)
    p.adjust.method <- match.arg(p.adjust.method)
#    n <- nrow(y)
    r <- t(apply(y, 1L, rank))
    q <- rank(apply(y, 1, function(u) max(u) - min(u)))
    s <- q * (r - (k+1)/2)
    w <- q * r
    ## s is a matrix of ranks within blocks (minus the average rank)
    ## multiplied by the ranked ranges of the blocks
    A <- sum(s^2)
    B <- sum(colSums(s)^2) / b
    S <- colSums(s)
    W <- colSums(w)
    if (dist == "TDist") {
        METHOD <- "posthoc-Quade test with TDist approximation"
        denom <- sqrt((2 * b * (A - B))/
                      ((b-1) * (k-1)))
        compare.stats <- function(i,j) {
            dif <- abs(S[i] - S[j]) 
            tval <- dif / denom
            return(tval)
        }
        PSTAT <- pairwise.table(compare.stats,levels(groups),
                                p.adjust.method="none")
        compare.levels <- function(i,j) {
            dif <- abs(S[i] - S[j]) 
            tval <- dif / denom
            pval <- pval <- 2 * pt(abs(tval),
                                   df=(b-1)*(k-1),
                                   lower.tail=FALSE)
            return(pval)
        }
        PVAL <- pairwise.table(compare.levels,levels(groups),
                               p.adjust.method=p.adjust.method)
    } else {
         METHOD <- paste("posthoc-Quade test",
                         "with standard-normal approximation",
                         sep="\t")
         n <- b * k
         denom <- sqrt((k * (k + 1) * (2 * n + 1) * (k-1))/
                           (18 * n * (n + 1)))
         nn <- length(w[,1])
         ff <- 1 / (nn * (nn + 1)/2)
         compare.stats <- function(i,j) {
            dif <- abs(W[i] * ff - W[j] * ff) 
            zval <- dif / denom
            return(zval)
        }
        PSTAT <- pairwise.table(compare.stats,levels(groups),
                                p.adjust.method="none")
        compare.levels <- function(i,j) {
            dif <- abs(W[i] * ff - W[j] * ff) 
            zval <- dif / denom
            pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
            return(pval)
        }
        PVAL <- pairwise.table(compare.levels,levels(groups),
                               p.adjust.method=p.adjust.method)
    }
    colnames(PSTAT) <- GRPNAMES[1:(k-1)]
    rownames(PSTAT) <- GRPNAMES[2:k]
    colnames(PVAL) <- GRPNAMES[1:(k-1)]
    rownames(PVAL) <- GRPNAMES[2:k]
    ans <- list(statistic = PSTAT,
                   p.value = PVAL,
                   method = METHOD,
                   p.adjust.method = p.adjust.method,
                   data.name = DNAME)
    class(ans) <- "PMCMR"
    ans
}
