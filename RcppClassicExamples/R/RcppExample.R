## RcppExample.R: Rcpp R/C++ interface class library example
##
## Copyright (C) 2008        Dirk Eddelbuettel
## Copyright (C) 2009 - 2012 Dirk Eddelbuettel and Romain Francois
##
## This file is part of Rcpp.
##
## Rcpp is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## Rcpp is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

RcppExample <- function(params, nlist, numvec, nummat, df, datevec, stringvec,
                        fnvec, fnlist) {

    ## Most of the input parameter checking here is not really
    ## necessary because it is done in the Rcpp code.

    ## Check that params is properly formatted.
    if (missing(params)) {
        cat("Setting default argument for params\n")
        params <- list(method='BFGS',
                       tolerance=1.0e-8,
                       maxIter=1000,
                       startDate=as.Date('2006-7-15'))
    }

    ## Check nlist
    if (missing(nlist)) {
        cat("Setting default argument for nlist\n")
        nlist <- list(ibm = 80.50, hp = 53.64, c = 45.41)
    } else if (!is.numeric(unlist(nlist))) {
        stop("The values in nlist must be numeric")
    }

    ## Check numvec
    if (missing(numvec)) {
        cat("Setting default argument for numvec\n")
        numvec <- seq(1,5) 			# numerical vector
    } else if (!is.vector(numvec)) {         ## Check numvec argument
        stop("numvec must be a vector");
    }

    ## Check nummat
    if (missing(nummat)) {
        cat("Setting default argument for nummat\n")
        nummat <- matrix(seq(1,20),4,5) # numerical matrix
    } else if (!is.matrix(nummat)) {
        stop("nummat must be a matrix");
    }

    ## Check df
    if (missing(df)) {
        cat("Setting default argument for data frame\n")
        df <- data.frame(a=c(TRUE, TRUE, FALSE), b=I(c('a','b','c')))
    }

    ## Check datevec
    if (missing(datevec)) {
        cat("Setting default argument for date vector\n")
        datestr <- c('2006-6-10', '2006-7-12', '2006-8-10')
        datevec <- as.Date(datestr, "%Y-%m-%d") # date vector
    }

    ## Check stringvec
    if (missing(stringvec)) {
        cat("Setting default argument for string vector\n")
        stringvec <- c("hello", "world", "fractal") # string vector
    }

    ## Check fnvec
    if (missing(fnvec)) {
        cat("Setting default argument for function vector\n")
        fnvec <- function(x) { sum(x) } # Add up components of vector
    }

    ## Check fnlist
    if (missing(fnlist)) {
        cat("Setting default argument for function list\n")
        fnlist <- function(l) { # Return vector with 1 added to each component
            vec <- c(l$alpha + 1, l$beta + 1, l$gamma + 1)
            vec
        }
    }

    ## Finally ready to make the call...
    val <- .Call("Rcpp_Example", params, nlist, numvec, nummat,
                 df, datevec, stringvec, fnvec, fnlist,
                 PACKAGE="RcppClassicExamples"
                 )

    ## Define a class for the return value so we can control what gets
    ## printed when a variable assigned this value is typed on a line by itself.
    ## This has the effect of calling the function print.RcppExample(). The
    ## function (defined below) simply prints the names of the fields that are
    ## available. Access each field with val$name.
    class(val) <- "RcppExample"

    val
}

print.RcppExample <- function(x,...) {
    cat('\nIn R, names defined in RcppExample return list:\n')
    cat('(Use result$name or result[[i]] to access)\n')
    namevec <- names(x)
    for(i in 1:length(namevec)) {
        cat(format(i, width=2), ': ', format(namevec[i], width=12))
        if (is.atomic(x[[i]])) {
            cat(format(x[[i]]))
        } else {
            cat(format("..."))
        }
        cat('\n')
    }
}

