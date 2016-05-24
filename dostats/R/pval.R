{###############################################################################
# pval.R
# Copyright 2012 Andrew Redd
# This file is part of the R package dostats
# Date: 5/30/2012
# 
# DESCRIPTION
# ===========
# Convenience functions for extracting pvalues and formatting test results.
# 
# LICENSE
# ========
# dostats is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# dostats is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# dostats. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################

#' Extract a p-value fomr a test result.
#' 
#' @param x a testing result object
#' @param extended should an extended result be given or a single p-value.
#' @param ... extra agruments
#' 
#' @details
#' This is a generic helper function for extracting p values from objects.
#' The idea being to extract the overall p-value for the model that can
#' be interpreted simply. 
#' 
#' 
#' @return either a single value (\code{extended=FALSE}) representing the 
#'         p-value of the test or a single row.
#'         \link{data.frame} object that also incldues extra information such as
#'         
#' @export
pval <- function(x, extended=F, ...)UseMethod('pval')

#'@export
pval.htest <- function(x, extended=F, ...) {
    if(!extended) {
        return(x$p.value)
    } else {
        with(x, data.frame(t=statistic, df=parameter, p.value))
    }
}

#'@export
pval.lm <- function(x, extended=F,...){
    tab <- anova(x)
    if(extended){
        as.data.frame(tab[1,])
    } else {
        tab[1, 5]
    }
}

#'@export
pval.default <- function(x, extended=F, ...){
    if('p.value' %in% names(x)) return(x$p.value)
    if(.hasSlot(x, 'p.value')) return(x@p.value)
    stop(sprintf("pval does not know how to extract a p-value from an object of class '%s'"
                , first(class(x))))
}


#' @export
#' @importFrom stats t.test
t.test.data.frame <- function(x, ...){
    message('Using dostats altered t.test for data.frames, you might consider using lm.')
    t.test(x[[1]] ~ x[[2]])
}




