## ################################################################################
## Program : NPST
## Language: R
## This program generalizes Hewitt's and Rogerson's
## nonparametric seasonality tests. Please see documentation for
## further details.
##
## Copyright (C) 2004-2014  Roland Rau
## 
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## ################################################################################

npst <- function(indata=NULL, long=12, peak=6, repts=100000, whole.distribution=FALSE,
                 siglevels=c(0.001,0.01,0.05,0.1), PARALLEL=FALSE, nodes=1) {
  # Function to generalize the nonparametric seasonality as provided
  # by Hewitt et al. (1971) and its (first) generalization by Rogerson
  # (1996)

  ################################################################################
  ################################################################################
  # Input:
  ################################################################################
  # indata              A numeric vector whose elements are the
  #  (numeric vector)   empirical number of events (e.g. deaths). The
  #                     length of the data is typically 12 (=months),
  #                     52 or 53 (weeks), or 365 or 366 (days).
  #                     Not providing 'indata'  is also okay (slightly
  #                     different output then).
  # 
  # long                The basic length of the data analyzed, i.e. if
  #  (numeric scalar)   we have monthly data, it would be 12. If
  #                     'indata' are provided, argument 'long' is
  #                     calculated based on argument 'indata'
  #
  # peak                Length of peak period. For instance, if we
  #  (numeric scalar)   assume that the 'peak season' lasts six months
  #                     for monthly data, peak would be six (=default value)
  #
  # repts               How many Monte Carlo simulation runs should be
  # (numeric scalar)    conducted?
  #
  # whole.distribution  Argument 'whole.distribution' indicates whether 
  #  (logical scalar)   the whole distribution should be returned (=TRUE)
  #                     or only the critical values (=FALSE).
  #
  # siglevels           For which significance levels should the
  #  (numeric vector)   corresponding required rank sums be
  #                     returned. Default settings are the 'typical' 
  #                     signficance levels of 0.001,0.01, 0.05, and
  #                     0.1
  # PARALLEL            Specify whether parallel estimation should be
  #  (boolean scalar)   performed (=TRUE) or not (FALSE=default).
  #
  # nodes               Specify on how many nodes the estimation
  #  (integer scalar)   should run (default=1). Only active if
  #                     argument PARALLEL=TRUE.
  #  
  ################################################################################
  ################################################################################

  # Output
  ################################################################################
  # maximum.rank.sum     The maximum rank sum theoretically possible
  # (numeric scalar)     with the given data.
  #
  # observed             The observed maximum rank sum (with the given data).
  # (numeric scalar)
  # 
  # observed.p.value     What is the p-value corresponding to the
  #  (numeric scalar)    observed maximum rank sum?
  #
  # critical.values      What are the required rank sums for the
  # (numeric data.frame) entered significance levels?
  #
  #        ########### ONLY IF whole.distribution=TRUE ###########
  # distribution         A data.frame specifying all possible
  # (numeric data.frame) rank-sums and their associated p-values.
  ################################################################################
  ################################################################################
  ################################################################################
  ################################################################################


  ### npst requires packages for parallel computation which only exist
  ### on Unix-like operating systems. If another operating system is
  ### used (e.g. MS Windows) and the argument PARALLEL is set to TRUE,
  ### it is overwritten to FALSE and a warning is issued.
  if ((.Platform$OS.type != "unix") && (PARALLEL==TRUE)) {
    PARALLEL <- FALSE
    cat("WARNING: Parallel execution not possible on your computing platform\n")
    cat("with package 'npst'. It requires packages which are currently only\n")
    cat("available on unix-like operating system.\n")
    cat("Argument PARALLEL set to FALSE.\n")
    cat("Execution of code being conducted on single core...\n\n")
  }
  
  ### Checking of Input:

  if (!is.null(indata)) {
    if (!is.vector(indata) || !is.numeric(indata)) {
      stop("Argument 'indata' needs to be a numeric vector!")
    } else {
      long <- length(indata)
    }
  }
  if ( (!is.numeric(long)) || (length(long)!=1)) {
    stop("Argument 'long' needs to be a numeric scalar!")
  }
  if ( (!is.numeric(peak)) || (length(peak)!=1)) {
    stop("Argument 'peak' needs to be a numeric scalar!")
  }

  if ( peak >= long) {
    stop("The value of 'long' needs to larger than the value of 'peak'!")
  }

  if ( (!is.numeric(repts)) || (length(repts)!=1)) {
    stop("Argument 'repts' needs to be a numeric scalar!")
  }

  if ( (length(whole.distribution) != 1) || (!is.logical(whole.distribution))) {
    stop("Argument 'whole.distribution' needs to be a logical scalar!")
  }
  if (!is.vector(siglevels) || !is.numeric(siglevels)) {
    stop("Argument 'siglevels' needs to be a numeric vector!")
  }
  ### END OF INPUT CHECKS

  # BEGIN of Calculating the observed maximum rank sum
  emp.rank <- function(indata, peak, returnvalue=0) {
    if (length(indata) > peak) {
      current.sum <- sum(indata[1:peak])
      if (current.sum > returnvalue) {
        returnvalue <- current.sum
        return(emp.rank(indata=indata[-1], peak=peak, returnvalue=returnvalue))
      } else {
        return(emp.rank(indata=indata[-1], peak=peak, returnvalue=returnvalue))
      }
    } else {
      return(returnvalue)
    }
  }
  emp.rank.wrapper <- function(indata, peak) {
    ranked.data <- rep(rank(indata),2)
    return(emp.rank(indata=ranked.data, peak=peak))
  }
  
  if (is.null(indata)) {
    actual.rank.sum <-
      "No Empirical Data were Provided => No Empirical Rank Sum Can Be Calculated"
  } else {
    actual.rank.sum <- emp.rank.wrapper(indata=indata, peak=peak)
  }
  # END of Calculating the observed maximum rank sum
  
  
  # what is the THEORETICAL maximum rank sum, given the input da?
  maxranksum <- sum(long:(long-peak+1))
  
  # initialize vector where results are stored with zeros.
  resvec <- rep(0, maxranksum)
  
  # initialize base data and repeat the values.
  basedata <- 1:long
  basedata2 <- rep(1:long,2)
  
  # calling external Fortran (90/95) routine to do the heavy work.
  if (PARALLEL) {### if boolean switch is TRUE, perform Monte Carlo
                 ### Simulation on as many cores as possible.
    require(parallel)
    no.of.workers <- nodes
    cat(paste(no.of.workers, "Nodes Will Be Used\n"))
    repts <- ceiling(repts/no.of.workers)
    
    results <- mclapply(X=1:no.of.workers,
                        FUN=function(x) .Fortran("npst", PACKAGE="npst",
                          lng=as.integer(long),
                          pk=as.integer(peak),
                          rpts=as.integer(repts),
                          mxrksum=as.integer(maxranksum),
                          resvec=as.integer(resvec),
                          basedata=as.integer(basedata),
                          basedata2=as.integer(basedata2))$resvec, mc.cores=nodes)
    results2 <- matrix(unlist(results), nrow=maxranksum,
                       ncol=no.of.workers, byrow=FALSE)
    pvals <- cumsum(apply(X=results2, MARGIN=1, FUN=sum)[order(maxranksum:1)]) / sum(results2)
  } else { ### from now on the previous code. In case when PARALLEL is
           ### not true, Monte Carlo Simulation will be performed only
           ### on one core

    # calling external Fortran (90/95) routine to do the heavy work.
    outsourcedloop <- .Fortran("npst", PACKAGE="npst",
                               lng=as.integer(long),
                               pk=as.integer(peak),
                               rpts=as.integer(repts),
                               mxrksum=as.integer(maxranksum),
                               resvec=as.integer(resvec),
                               basedata=as.integer(basedata),
                               basedata2=as.integer(basedata2))

    # Calculate the p-values, i.e. the cumulative distribution function.
    pvals <- cumsum(outsourcedloop$resvec[order(maxranksum:1)])/repts
  }
  
  # Calculate the p-values, i.e. the cumulative distribution function.
  maxisums <- maxranksum:1
  returndf <- data.frame(Max.Rank.Sum=maxisums, p.values=pvals)
  rownames(returndf) <- 1:(nrow(returndf))

  # Calculating the positions (=which rows) in the data.frame
  # 'returndf' where the significance levels provided by the user are
  # reached.
  the.positions <- NULL
  for (i in 1:(length(siglevels))) {
    current.pos <- which.min( abs(pvals - siglevels[i]))
    current.posis <- (current.pos-2) : (current.pos+2)
    the.positions <- c(the.positions,current.posis)
  }
  the.positions <- sort(unique(the.positions[(the.positions > 0) & (the.positions <= (length(pvals)))]))
  

  ## Adjusting the final output -- depending on the input parameters
  if (whole.distribution) {
    if (is.null(indata)) {
      returnlist <- list(maximum.rank.sum=maxranksum,
                         critical.values=returndf[the.positions,],
                         distribution=returndf)
      
    }
    if (!is.null(indata)) {
      returnlist <- list(maximum.rank.sum=maxranksum,
                         observed=actual.rank.sum,
                         observed.p.value=pvals[maxranksum-actual.rank.sum+1],
                         critical.values=returndf[the.positions,],
                         distribution=returndf)
    }
  }
  if (!whole.distribution) {
    if (is.null(indata)) {
      returnlist <- list(maximum.rank.sum=maxranksum,
                         critical.values=returndf[the.positions,])
    }
    if (!is.null(indata)) {    
      returnlist <- list(maximum.rank.sum=maxranksum,
                         observed=actual.rank.sum,
                         observed.p.value=pvals[maxranksum-actual.rank.sum+1],
                         critical.values=returndf[the.positions,])
    }
  }
  return(returnlist)
}
