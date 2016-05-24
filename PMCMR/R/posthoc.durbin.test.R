#  posthoc.durbin.test.R
#  Part of the R package PMCMR
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

posthoc.durbin.test <- function(y, ...) UseMethod("posthoc.durbin.test")

posthoc.durbin.test.default <-
function(y, groups, blocks,  p.adjust.method = p.adjust.methods, ...)
{
    DNAME <- deparse(substitute(y))

    if (is.matrix(y)) {
        GRPNAME <- colnames(y)
        g <- factor(c(col(y)))
        b <- factor(c(row(y)))
        yy <- as.vector(y)
        datf1 <- data.frame(yy, b, g)
        datf2 <- datf1[!is.na(datf1[,1]),]
        blocks <- factor(datf2[,2])
        groups <- factor(datf2[,3])
        y <- datf2[,1]
        ## Clean up
        rm(b, g, yy, datf1, datf2)
    }
    else {
        if (anyNA(groups) || anyNA(blocks))
            stop("NA's are not allowed in 'groups' or 'blocks'")
        if (any(diff(c(length(y), length(groups), length(blocks))) != 0L))
            stop("'y', 'groups' and 'blocks' must have the same length")
        DNAME <- paste(DNAME, ", ", deparse(substitute(groups)),
                       " and ", deparse(substitute(blocks)), sep = "")
        groups <- factor(groups)
        blocks <- factor(blocks)
        GRPNAME <- levels(groups)
    }
    
    ## Need to ensure consistent order.
    o <- order(blocks, groups)
    y <- y[o]
    groups <- groups[o]
    blocks <- blocks[o]
    
    p.adjust.method = match.arg(p.adjust.method)
    t <- nlevels(groups)
    b <- nlevels(blocks)
    r <- unique(table(groups))
    k <- unique(table(blocks)) 
    rij <- unlist(tapply(y, blocks, rank))
    Rj <- tapply(rij, groups, sum)
    ## Taken from NIST
    A <- sum(rij^2)
    C <- (b * k * (k + 1)^2) / 4
    D <- sum(Rj^2) - r * C
    T1 <- (t - 1) / (A - C) * D

    denom <- sqrt(((A - C) * 2 * r) / (b * k - b - t + 1) *
                      (1 - T1 / (b * (k -1))))
    df <- b * k - b - t + 1
   # Pairwise comparisons
    compare.stats <- function(i,j) {
        dif <- abs(Rj[i] - Rj[j]) 
        tval <- dif / denom
        return(tval)
    }
    PSTAT <- pairwise.table(compare.stats,levels(groups),
                            p.adjust.method="none" )
    
    compare.levels <- function(i,j) {
        dif <- abs(Rj[i] - Rj[j]) 
        tval <- dif / denom
        pval <- 2 * pt(q=abs(tval), df=df, lower.tail=FALSE)
        return(pval)
    }
    PVAL <- pairwise.table(compare.levels,levels(groups),
                           p.adjust.method=p.adjust.method)

    METHOD <- paste("Durbin's test for a two-way", "
                 balanced incomplete block design", sep="\t")
    colnames(PSTAT) <- GRPNAME[1:(t-1)]
    rownames(PSTAT) <- GRPNAME[2:t]
    colnames(PVAL) <- colnames(PSTAT)
    rownames(PVAL) <- rownames(PSTAT)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method =p.adjust.method)
    class(ans) <- "PMCMR"
    ans
}
