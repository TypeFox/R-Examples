### ----------------------------------------------------------------------------
### This file is part of boolean3
###
### Copyright (C) 2011--2014 Jason W. Morgan <morgan.746@osu.edu>
###
### boolean3 represents a substantial re-write of the original boolean package
### developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
### was developed under the direction of Bear Braumoeller and with support from
### The Ohio State University's College of Social and Behavioral Sciences.
###
### boolean3 and is free software: you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the Free
### Software Foundation, either version 3 of the License, or (at your option)
### any later version.
###
### This program is distributed in the hope that it will be useful, but WITHOUT
### ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
### FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
### more details.
###
### You should have received a copy of the GNU General Public License along with
### this program.  If not, see <http://www.gnu.org/licenses/>.
###
### ----------------------------------------------------------------------------


## =============================================================================
## NOTE: I no longer think this is a good way to deal with WARNINGS and
## ERRORS. These should be incorporated directly into the code to facilitate
## debugging.
## =============================================================================

## Prints \code{warning} message given specified code.
##
## This function prints a \code{warning} message given a specified code. If no
## such code exists, it returns a generic message.
## 
## @param msg message code, should be an integer. 
## @return prints a \code{warning}.
## @author Jason W. Morgan <morgan.746@@osu.edu>
warn_msg <- function(msg) {
  warn.msg <- switch(as.character(msg),
                     ## boolboot warnings.
                     "61" = "convergence error in bootstrap (code 61)",
                     "62" = "number of bootstrap samples is probably too small",
                     "63" = "number of errors in convergence met exceeded 10% of the \n  number of bootstrap samples. Estimate likely unstable."
                     )

  if(is.null(warn.msg))
    warn.msg <- paste("Undefined warning (", msg, ")", sep = "")

  warning(warn.msg, call. = FALSE)
}

## Prints convergence \code{warning} message given specified code.
##
## This function prints a \code{warning} message given a specified convergence
## code error. If no such code exists, it returns a generic message.
##
## @param id message id.
## @param msg detailed error message. Not printed if NULL. 
## @return prints a \code{warning}.
## @author Jason W. Morgan <morgan.746@@osu.edu>
conv_msg <- function(id, msg=NULL, method=NULL) {
  warn.msg <- switch(as.character(id, msg),
                     ## optimx generic error messages
                     "1" = "iteration limit 'maxit' had been reached (optimx code 1)",
                     "20" = "initial parameter set is inadmissible (optimx code 20)",
                     "21" = "intermediate parameters set is inadmissible (optimx code 21)",
                     "10" = "degeneracy of the Nelder-Mead simplex (optimx code 10)",
                     "51" = "warning from 'L-BFGS-B' method (optimx code 51)",
                     "52" = "error from the 'L-BFGS-B' method (optimx code 52)"
                     )

  if(is.null(warn.msg))
    warn.msg <- paste("Undefined convergence warning (", msg, ")", sep = "")

  warning(paste(method, ": ", warn.msg, sep=""), call.=FALSE, immediate.=FALSE)
  if (!is.null(msg))
    warning(paste(method, ": ", msg, sep=""), call.=FALSE, immediate.=FALSE)
}

## Prints an error message and calls \code{stop}.
##
## This function prints an error message given a specific code and then calls
## \code{stop}.
## 
## @param msg message code, should be an integer. 
## @return halts execution.
## @author Jason W. Morgan <morgan.746@@osu.edu>
error_msg <- function(msg) {
  err.msg <- switch(as.character(msg),
                    ## boolprof, boolprob, and boolboot error codes.
                    "22" = "boolboot requires a boolean object as outputed by boolprep (code 22)",
                    "23" = "model object already contains bootstrapping results. Halting to prevent data loss (code 23)",
                    )
  
  if(is.null(err.msg))
    err.msg <- paste("Undefined fatal error (", msg, ")", sep = "")
  stop(err.msg, call.=FALSE)
}

misc_msg <- function(msg) {
  msg <- switch(as.character(msg),
                ## boolboot messages
                "1"  = "Warning: boolean likelihood functions are often *highly* irregular,\n   thus normal-theory p-values are suspect. Bootstrapped confidence\n   intervals should be used for inference.\n",
                "61" = "Warning: number of bootstrap samples is probably too small. Try \n         increasing 'n' in the call to boolboot.",
                "62" = "Warning: number of errors in convergence exceeded 10% of the \n         number of bootstrap samples. Estimates likely unstable."
                )
  
  if(is.null(msg))
    msg <- paste("Undefined message (", msg, ")", sep = "")
  message(msg)
}
