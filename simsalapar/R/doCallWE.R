## Copyright (C) 2012 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 2 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


mkTimer <- function(gcFirst) {
    if(gcFirst)
	 function(expr) 1000 * system.time(expr, gcFirst=TRUE)[[1]]
    else function(expr) 1000 * system.time(expr, gcFirst=FALSE)[[1]]
}

##' @title Safety Wrapper for Innermost Computations
##' @param f a function, given data + parameters, computes statistic
##' @param argl list of arguments for f()
##' @param timer a function like \code{\link{system.time}()}, by default,
##'              measure user time in milliseconds
##' @return list with components:
##'	    value:   f(<argl>) if that worked, NULL otherwise
##'	    error:   error message (of class simpleError) or NULL
##'	    warning: warning message (of class simpleWarning) or NULL
##'	    time:    time measured by timer()
##' @author Marius Hofert, Martin Maechler
##' { doCallWE }
doCallWE <- function(f, argl, timer = mkTimer(gcFirst=FALSE))
{
    tim <- timer( res <- tryCatch.W.E( do.call(f, argl) )) # compute f(<argl>)
    is.err <- is(val <- res$value, "simpleError") # logical indicating an error
    list(value   = if(is.err) NULL else val, # value (or NULL in case of error)
	 error   = if(is.err) val else NULL, # error (or NULL if okay)
	 warning = res$warning, # warning (or NULL)
	 time    = tim) # time
}
##' { end }
