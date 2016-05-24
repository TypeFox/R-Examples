{###############################################################################
# wargs.R
# This file is part of the R package dostats.
# 
# Copyright 2012 Andrew Redd
# Date: 6/1/2012
# 
# DESCRIPTION
# ===========
# unit tests for recombine function
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
wrap_function <- function(symb, args, envir, attrs=NULL){
#' @param symb symbolic name of function
#' @param args pairlist of arguments to use
#' @param envir the environment for the function
#' @param attrs attributes to carry forward
    newargs <- args
    for(a in setdiff(names(newargs), '...')) {
        newargs[[a]] <- as.name(a)
    }
    wdots <- names(newargs)=='...'
    newargs[wdots] <- list(as.symbol('...'))
    names(newargs)[wdots] <- ''
    c1 <- as.call(append(list(symb), as.list(newargs)))
    c2 <- as.call(c(as.name('{'), (c1)))
    newfunction <- as.function(append(args, list(c2)), envir=envir)
    if(!is.null(attrs))
        attributes(newfunction) <- attrs
    return(newfunction)
}

#'  Call with arguments.
#' 
#'  @param f a function
#'  @param ... extra arguments
#'  @param args alternate way to provide arguments as a pairlist.
#'  @param envir environment to use for the function.
#' 
#'  @return a function that takes 1 argument and calls f with the 
#'  single argument and the additional \code{...} appended.
#'  @export
#'  @keywords utilities, misc
#'  @examples 
#'  mean2 <- wargs(mean, na.rm=TRUE)
wargs <- function(f, ..., args=pairlist(...), envir = parent.frame()){
    symb <- substitute(f)
    af   <- formals(base::args(f))
    new.args <- c(af[setdiff(names(af), names(args))], args)
    wrap_function(symb, new.args, envir, attrs = attributes(f))
}



#' Create a function that redirects to the named function.
#' 
#' This is usefull for debugging to know what function has been called 
#' form within do.call or plyr functions.
#' 
#' @param f     a function to wrap a call around
#' @param envir environment to use for the function.
#'
#' @export
redirf <- function(f, envir=parent.frame()){
    symb <- substitute(f)
    args <- formals(f)
    wrap_function(symb, args, envir)
}
