## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015 - 2016  Stefan Engblom
## Copyright (C) 2015 - 2016  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' Check integer arguments
##'
##' Raise an error if any of the arguments are non-integer.
##' @param len Expected length of the infectious pressure vector
##' @param ... The arguments to check
##' @keywords internal
##' @return invisible(NULL)
check_infectious_pressure_arg <- function(len, ...) {
    arg <- list(...)
    for (i in seq_len(length(arg))) {
        if (!is.numeric(arg[[i]])) {
            stop(paste0("Invalid '",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "': must be numeric vector"))
        }

        if (!is.null(dim(arg[[i]]))) {
            stop(paste0("Invalid '",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "': must be numeric vector"))
        }

        if (!identical(length(arg[[i]]), len)) {
            stop(paste0("Invalid '",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "': must be numeric vector with length 'nrow(u0)'"))
        }

        if (any(arg[[i]] < 0)) {
            stop(paste0("Invalid '",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "': must be numeric vector with non-negative values"))
        }
    }

    invisible(NULL)
}

##' Check integer arguments
##'
##' Raise an error if any of the arguments are non-integer.
##' @param ... The arguments to check
##' @keywords internal
##' @return invisible(NULL)
check_integer_arg <- function(...) {
    arg <- list(...)
    for (i in seq_len(length(arg))) {
        if (is.null(arg[[i]])) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "' is missing"))
        }

        if (!is.numeric(arg[[i]])) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "' must be integer"))
        }

        if (!all(is_wholenumber(arg[[i]]))) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "' must be integer"))
        }
    }

    invisible(NULL)
}

##' Check arguments for 'gdata'
##'
##' Raise an error if any of the arguments are not ok.
##' @param ... The arguments to check
##' @keywords internal
##' @return invisible(NULL)
check_gdata_arg <- function(...) {
    arg <- list(...)
    for (i in seq_len(length(arg))) {
        if (is.null(arg[[i]])) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "' is missing"))
        }

        if (!is.numeric(arg[[i]])) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "' must be numeric"))
        }

        if (!identical(length(arg[[i]]), 1L)) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "' must be of length 1"))
        }
    }

    invisible(NULL)
}

##' Check arguments for interval endpoints
##'
##' Raise an error if any of the arguments are not ok.
##' @param len Exprected length of each of the interval endpoint
##' vectors
##' @param ... The arguments to check
##' @keywords internal
##' @return invisible(NULL)
check_end_t_arg <- function(len, ...) {
    arg <- list(...)
    names(arg) <- match.call(expand.dots = FALSE)$'...'

    for (i in seq_len(length(arg))) {
        if (!identical(length(arg[[i]]), len)) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "' must be of length 1 or 'nrow(u0)'"))
        }
    }

    ## Check interval endpoints
    if (!all(0 <= arg$end_t1))
        stop("'end_t1' must be greater than or equal to '0'")
    if (!all(arg$end_t1 < arg$end_t2))
        stop("'end_t1' must be less than 'end_t2'")
    if (!all(arg$end_t2 < arg$end_t3))
        stop("'end_t2' must be less than 'end_t3'")
    if (!all(arg$end_t3 < 364))
        stop("'end_t3' must be less than '364'")
    if (!all(0 <= arg$end_t4))
        stop("'end_t4' must be greater than or equal to '0'")
    if (!all(arg$end_t4 <= 365))
        stop("'end_t4' must be less than or equal to '365'")
    if (!all((arg$end_t4 < arg$end_t1) | (arg$end_t3 < arg$end_t4)))
        stop("'end_t4' must be less than 'end_t1' or greater than 'end_t3'")

    invisible(NULL)
}
