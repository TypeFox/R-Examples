#  durbin.test.R
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

durbin.test <- function(y, ...) UseMethod("durbin.test")

durbin.test.default <-
function(y, groups, blocks, ...)
{
    DNAME <- deparse(substitute(y))

    if (is.matrix(y)) {
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
    }
    
    ## Need to ensure consistent order.
    o <- order(blocks, groups)
    y <- y[o]
    groups <- groups[o]
    blocks <- blocks[o]

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
    
    STATISTIC <- T1
    PARAMETER <- t - 1
    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
    names(STATISTIC) <- "Durbin chi-squared"
    names(PARAMETER) <- "df"

    structure(list(statistic = STATISTIC,
                   parameter = PARAMETER,
                   p.value = PVAL,
                   method = "Durbin rank sum test",
                   data.name = DNAME),
              class = "htest")
}

durbin.test.formula <-
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
    y <- do.call("durbin.test", as.list(mf))
    y$data.name <- DNAME
    y
}
