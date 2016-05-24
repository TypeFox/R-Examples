#*******************************************************************************
#
# Local Approximate Gaussian Process Regression
# Copyright (C) 2013, The University of Chicago
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
#
#*******************************************************************************
 

## getDs
##
## calculate an initial starting value and maximum value

getDs <- function(X, p=0.1, samp.size=1000)
  {
    if(nrow(X) > samp.size) {
      i <- sample(1:nrow(X), samp.size)
      X <- X[i,]
    }
    D <- distance(X)
    D <- D[upper.tri(D)]
    D <- D[D > 0]
    dstart <- as.numeric(quantile(D, p=p))
    ## diag(D) <- NA
    ## dstart <- (mean(apply(D, 1, min, na.rm=TRUE)) + mean(D, na.rm=TRUE))/2
    dmax <- max(D, na.rm=TRUE)
    dmin <- min(D, na.rm=TRUE)
    return(list(start=dstart, min=dmin, max=dmax))
  }


## check.arg:
## 
## check a gamma-defined argument as created by garg and 
## darg beloe

check.arg <- function(d)
  {
    if(length(d$max) != 1 || d$max <= 0)
      stop("d$max should be a positive scalar")
    if(length(d$min) != 1 || d$min < 0 || d$min > d$max)
      stop("d$min should be a postive scalar < d$max")
    if(any(d$start < d$min) || any(d$start > d$max))
      stop("(all) starting d-value(s) should be a positive scalar in [d$min, d$max]")
    if(length(d$mle) != 1 || !is.logical(d$mle))
      stop("d$mle should be a scalar logical")
    if(length(d$ab) != 2 || any(d$ab < 0)) 
      stop("$ab should be a positive 2-vector")
    invisible(NULL)
  }


## garg:
##
## process the g argument to approxGP and localGP, which
## specifies starting values, ranges, whether or not mle calcuations
## should be made, and priors

garg <- function(g, y)
  {
    ## coerce inputs
    if(is.null(g)) g <- list()
    else if(is.numeric(g)) g <- list(start=g)
    if(!is.list(g)) stop("g should be a list or NULL")

    ## check mle
    if(is.null(g$mle)) g$mle <- FALSE
    if(length(g$mle) != 1 || !is.logical(g$mle))
      stop("g$mle should be a scalar logical")

    ## check for starting and max values
    if(is.null(g$start) || (g$mle && (is.null(g$max) || is.null(g$ab) || is.na(g$ab[2]))))
      r2s <- (y - mean(y))^2

    ## check for starting value  
    if(is.null(g$start)) g$start <- as.numeric(quantile(r2s, p=0.025))

    ## check for max value
    if(is.null(g$max)) {
      if(g$mle) g$max <- max(r2s)
      else g$max <- max(g$start)
    }

    ## check for dmin
    if(is.null(g$min)) g$min <- sqrt(.Machine$double.eps)

    ## check for priors
    if(!g$mle) g$ab <- c(0,0)
    else {
      if(is.null(g$ab)) g$ab <- c(3/2, NA)
      if(is.na(g$ab[2])) {
        s2max <- mean(r2s)
        g$ab[2] <- Igamma.inv(g$ab[1], 0.95*gamma(g$ab[1]), 1, 0)/s2max
      }
    }

    ## now check the validity of the values, and return
    check.arg(g)
    return(g)
  }


## darg:
#3
## process the d argument to approxGP and localGP, which
## specifies starting values, ranges, whether or not mle calcuations
## should be made, and priors

darg <- function(d, X, samp.size=1000)
  {
    ## coerce inputs
    if(is.null(d)) d <- list()
    else if(is.numeric(d)) d <- list(start=d)
    if(!is.list(d)) stop("d should be a list or NULL")

    ## check for mle
    if(is.null(d$mle)) d$mle <- TRUE

    ## check if we need to build Ds
    if(is.null(d$start) || (d$mle && 
      (is.null(d$max) || is.null(d$min) || is.null(d$ab) || is.na(d$ab[2]))))
      Ds <- getDs(X, samp.size=samp.size)

    ## check for starting value
    if(is.null(d$start)) d$start <- Ds$start

    ## check for max value
    if(is.null(d$max)) {
      if(d$mle) d$max <- Ds$max
      else d$max <- max(d$start)
    }

    ## check for dmin
    if(is.null(d$min)) {
      if(d$mle) d$min <- Ds$min/2
      else d$min <- min(d$start)
      if(d$min < sqrt(.Machine$double.eps))
        d$min <- sqrt(.Machine$double.eps)
    }

    ## check for priors
    if(!d$mle) d$ab <- c(0,0)
    else {
      if(is.null(d$ab)) d$ab <- c(3/2, NA)
      if(is.na(d$ab[2]))
        d$ab[2] <- Igamma.inv(d$ab[1], 0.95*gamma(d$ab[1]), 1, 0)/Ds$max
    }

    ## now check the validity of the values, and return
    check.arg(d)
    return(d)
  }



## Igamma.inv:
##
## calculate the beta parameter of an Inverse Gamma
## distribution with alpha parameter a at location
## y

Igamma.inv <- function(a, y, lower=FALSE, log=FALSE)
  {
    ## call the C routine
    r <- .C("Igamma_inv_R",
            a = as.double(a),
            y = as.double(y),
            lower = as.integer(lower),
            log = as.integer(log),
            result = double(1),
            PACKAGE = "laGP")

    return(r$result)
  }
