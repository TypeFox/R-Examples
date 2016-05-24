#******************************************************************************* 
#
# Particle Learning of Gaussian Processes
# Copyright (C) 2010, University of Cambridge
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
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## visible definition for R CMD CHECK
PL.env <- new.env()

## PL.clear:
##
## clear the pall and peach particle information

PL.clear <- function()
  {
    PL.env$pall <- list()
    PL.env$peach <- list()
    PL.env$psave <- NULL
  }


## initialize
PL.clear()


## resample:
##
## use the provided predprob function to do the PL-resample
## step by frist calculating the resample weights with
## lpredprob and then resampling the particles according
## to the weights

resample <- function(z, lpredprob, prior)
  {
    ## allocate space for the weights
    P <- length(PL.env$peach)
    lweights <- rep(NA, P)

    ## calculate the weight of each particle
    for(p in 1:P)
      lweights[p] <- lpredprob(z, PL.env$peach[[p]], prior)

    ## re-normalize the log weights
    weights <- renorm.lweights(lweights)
    
    ## resample particle indices according to the weighst
    indices <- sample(1:P, P, prob=weights, replace=TRUE)
    Pnew <- length(unique(indices))
    if(P > 1 && Pnew == 1) stop("total degeneration")
    ## cat("unique particles ", Pnew, "\n")

    ## re-allocate the list of particles according to
    ## the re-sampled indices
    newp <- list()
    for(p in 1:P) newp[[p]] <- PL.env$peach[[indices[p]]]

    ## copy in the new particles
    PL.env$peach <- newp
  }



## PL:
##
## the master particle learning function for via the propagation
## and resample steps--  a print statement can be made every
## progress particles (none if zero)

PL <- function(dstream, start, end, init, lpredprob, propagate, prior=NULL,
               addpall=NULL, params=NULL, save=NULL, P=100, progress=10,
               cont=FALSE, verb=1)
  {
    ## calculate the starting (initalization)
    ## if(start >= end) stop("no PL to do since start >= end")

    ## initialize if not picking up where we left off
    if(!cont) {
      
      ## clear the global variables
      PL.clear()
    
      ## push the starting data to the global
      if(!is.null(addpall)) addpall(dstream(1,start))
    
      ## do a few MH steps on the particles to get them to
      ## move from being samples from the prior to being samples from
      ## the posterior up to time start
      PL.env$peach <- list()
      PL.env$peach[[1]] <- init(NULL, prior=prior)
      if(P > 1) {
        for(i in 2:P) {
          PL.env$peach[[i]] <- init(PL.env$peach[[i-1]], prior=prior)
          ## print progress
          if(verb != 0) cat("  particle init: ", signif(100*i/P,2), "%   \r", sep="")
        }
      }
      ## cap off progress printing and possibly add a plot of params
      if(verb != 0) cat("\n")
      if(progress > 0 && !is.null(params)) {
        phist(start, params(), "init=")
        if(verb == -1) readline("init done, press RETURN to continue: ")
      }
    }

    ## no PL to do, only MCMC
    if(start == end) return(PL.env$peach)
    
    ## do the particle learning
    for(t in (start+1):end) {

      ## print progress
      if(verb != 0)
        cat("  particle learning: ", signif(100*t/end,2), "%   \r", sep="")
      
      ## re-sample step according to the next input-output pair
      z <- dstream(t)
      resample(z, lpredprob, prior)  ## causes a change in peach

      ## optionally plot histograms of the parameters
      if((progress > 0) && (t %% progress == 0) && !is.null(params) && verb == -1) {
        phist(t, params(), "res=")
        readline("\nresample plotted, press RETURN for propogate: ")
      }
      
      ## add the t-he t-th sample global
      if(!is.null(addpall)) addpall(z)
      
      ## propagate each particle
      for(p in 1:P) PL.env$peach[[p]] <- propagate(z, PL.env$peach[[p]], prior)
      
      ## optionally plot histograms of the params
      if((progress > 0) && (t %% progress == 0) && !is.null(params)) {
        if(verb == -1) { str <- "prop=" } else { str <- "t=" }
        phist(t, params(), str)
        if(verb == -1) readline("propagate plotted, press RETURN to continue: ")
      }

      ## optionally save information about every particle
      if(!is.null(save)) save()
    }

    ## cap off progress meter
    if(progress > 0) cat("\n")

    ## return the PL particles
    return(PL.env$peach)
  }


## phist:
##
## simple function that calculates the layout for
## plotting histograms of the parameters in all particles
## stored in the data frame p.  t is an iteration number
## for the main title which also include the parameter
## ranges, and number of unique versions of each parameter
## all pre-ceeded by an optional string str

phist <- function(t, p, str=NULL, cpar=TRUE)
  {
    if(is.null(p)) return()

    if(cpar) {
      rows <- floor(sqrt(length(p)))
      cols <- floor(length(p) / rows)
      while(rows * cols < length(p)) cols <- cols + 1
      par(mfrow=c(rows, cols), bty="n")
    }
    for(i in 1:ncol(p)) {
      rd <- range(p[,i])
      tit <- paste(str, t, " r=(", signif(rd[1],2), ",", signif(rd[2], 2), ")",
                 " u=", length(unique(p[,i])), sep="")
      hist(p[,i], main=tit, xlab=names(p)[i])
    }
  }


## papply:
##
## apply the function fun to each particle and collect the output
## into a list, a print statement can be made every progress
## particles (none if zero)

papply <- function(fun, verb=1, pre="", ...)
  {
    P <- length(PL.env$peach)

    ## init the return list
    Ztp <- list()

    ## collect the summary statistics from the particles
    for(p in 1:P) {

      ## optional progress printing
      if(verb) cat(" ", pre, " applying: ", signif(100*p/P,2), "%   \r", sep="")

      ## predict at the XX locations for particle p
      Ztp[[p]] <- fun(Zt=PL.env$peach[[p]], ...)
    }

    ## finish off verbosity printing
    if(verb > 0) cat("\n")

    ## return the list
    return(Ztp)
  }


## renorm.lweights:
##
## re-normalize a (un-logged) weight vector

renorm.weights <- function(weights)
  {
    return(weights/sum(weights))
  }


## renorm.lweights:
##
## for numerical stability, re-normalizes a vector
## of log-weights by adding in a factor to make the
## (log) weights smaller, than returns a normalized
## version of the weights on the un-logged scale

renorm.lweights <- function(lweights)
  {
    lweights= lweights - max(lweights) 
    return(renorm.weights(exp(lweights)))

}

