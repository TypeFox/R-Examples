{###############################################################################
# utils.R
# Copyright 2012 Andrew Redd
# Date: 5/30/2012
# 
# DESCRIPTION
# ===========
# Convenience functions and utilities for use in dostats.
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


#' @rdname head-tail-utils.Rd
#' @name first
#' @title Head/Tail shortcuts
#' 
#' @description
#' Shortcuts for \code{head(x,1)} and \code{tail(x, 1)}
#' 
#' @param x   vector object
#' @param ... passed on to head or tail
#' @param n   the new number to take of only one.
#' 
#' @export
first <- wargs(head, n=1)
#' @rdname head-tail-utils.Rd
#' @export
last <- wargs(tail, n=1)

is_uniform <- function(x){
    all(sapply(x, isTRUE %.% all.equal, x[[1]]))
}

#' List rows of a data frame in a list.
#' 
#' @param d a data.frame
#' 
#' @export
listrows <- function(d){
    FUN <- data.frame
    dots <- as.list(d)
    MoreArgs <- list(stringsAsFactors=FALSE)
    .mapply(data.frame, dots, MoreArgs)
}

#' Make a helper ID counter
#' 
#' @param startat where to start counting
#' 
#' @export
make_new_id <- function(startat=0){
    ..next_id <- startat
    list(
    new = function(n=1){
        ..next_id <<- ..next_id +n
        return((..next_id-n+1):..next_id)
    },
    reset = function(at=0){..next_id <<- at},
    curr  = function(){..next_id}
    )
}

#' Make a call with extra arguments incorporated into call.
#'
#' Usefull for using with plyr functions
#'
#' @param args a list of arguments
#' @param ... extra arguments to be incorporated into args
#' @param what the function to execute
#' @param quote should the arguments be quoted
#' @param envir the environment to call the function in
#' 
#' @seealso \code{\link{do.call}} which this function wraps.
#' @export
make_call <- function(args, ..., what, quote = F, envir=parent.frame()){
    args <- append(args, list(...))
    do.call(what=what, args = args, quote=quote, envir=envir)
}

#' Return the current function
#' @seealso \code{\link{sys.function}}
#' @export
me<- function(){
    # print(sys.function(-2))
    i <- 0
    while(T){
        i <- i-1
        fun <- sys.function(i)
        if(identical(environment(fun), asNamespace("plyr"))) next
        if(identical(environment(fun), asNamespace("harvestr"))) next
        if(!identical(fun, me)) return(fun)
    }
}


#' Fill vector to length with a specified value
#' 
#' @param x vector
#' @param l length
#' @param with What to fill with
#' @param after where to insert
#' 
#' @export
fill_v <- function(x, l=length(x), with=last(x), after=length(x)){
    stopifnot(length(x) <= l)
    append(x, rep(with, l-length(x)), after=after)
}
