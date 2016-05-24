#******************************************************************************* 
#
# Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
# Copyright (C) 2005, University of California
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
# Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
#
#*******************************************************************************


## optim.ptgpf:
##
## find the minima of the MAP predictive (kriging) surface
## encoded in the tgp object, starting at the specified spot
## restrected to the provided rectangle with the specified
## method  -- eventually we should be calclating and using GP 
## derivative information

optim.ptgpf <- function(start, rect, tgp.obj,
                        method=c("L-BFGS-B", "Nelder-Mead", "BFGS",
                          "CG", "SANN", "optimize"))
{
  ## check the method argument
  method <- match.arg(method)

  ## ptgpf:
  ##
  ## predict at x for the MAP tgp object, to be used by optim
  ## for finding the minimum of the MAP kriging surface
  
  ptgpf <- function(x, tgp.obj, rect=NULL)
    {
      ## only need to check rectangle when 2-d or more
      if(!is.null(rect)) for(i in nrow(rect))
          if(x[i] < rect[i,1] || x[i] > rect[i,2]) return(Inf)
      
      ## necessary b/c check.matrix doesn't know correct ncol
      if(!is.null(rect)) x <- matrix(x, ncol=nrow(rect))
      
      ## run predict
      out <- predict(tgp.obj, XX=x, pred.n=FALSE)
      return(as.vector(out$ZZ.km))
    }
  
  ## optimize is for 1-d data only
  if(method == "optimize") {
    if(nrow(rect) != 1) ## check if optimize method is appropriate
      stop("method=\"optimize\" only valid for 1-d functions")
    opt <- optimize(ptgpf, interval=rect[1,], tgp.obj=tgp.obj)
    return(list(par=opt$minimum, value=opt$objective, convergence=1))
  }

  ## otherwise use optim in some way
  if(method == "L-BFGS-B") {  ## use the boundary informatoin in rect
    opt <- optim(par=start, ptgpf, method=method, tgp.obj=tgp.obj,
                 rect=rect, lower=rect[,1], upper=rect[,2])
    
  } else { ## otheriwise, apply a method without boundaries
    opt <- optim(par=start, ptgpf, method=method, tgp.obj=tgp.obj,
                 rect=rect)
  }
  
  ## return
  return(opt)
}


## tgp.cands:
##
## create NN candidate locations (XX) either via Latin Hypercube
## sample (LHS), or sequential treed D-optimal design (based on an
## initial LHS

tgp.cands <- function(tgp.obj, NN, cands=c("lhs", "tdopt"), rect, X, verb=0)
  {
    ## check the cands argument
    cands <- match.arg(cands)

    ## return a latin hypercibe sample
    if(cands == "lhs") return(lhs(NN, rect))

    ## return a sequential treed D-optimal sample from initial LHS cands
    Xcand <- lhs(10*NN, rect)
    if(is.null(tgp.obj)) XX <- dopt.gp(NN, X=X, Xcand, verb=verb)$XX
    else XX <- tgp.design(NN, Xcand, tgp.obj, verb=verb)
    XX <- matrix(XX, ncol=ncol(X))
    return(XX)
  }


## optim.tgp:

## execute one step in a search for the global optimum (minimum) of a
## noisy function (f) bounded in rect with starting (X,Z) data provided:
## fit a tgp model and predict creating NN+{1,2} candidates and
## select the one (or two) which give have the highest expected improv
## statistic.  NN of the candidates come from cands (lhs or tdopt),
## plus one which is the location of the minima found (e.g.,) via calling
## optim (with particular method) on the MAP btgpm predictive surface
## (passed in with prev).  When as != "none" an additional candidate
## is also selected, which has the highest expected alc or alm statistic
## The new X (which may be 1-3 rows) are returned

optim.step.tgp <-
  function(f, rect, model=btgp, prev=NULL, X=NULL, Z=NULL,
           NN=20*length(rect), improv=c(1,5), cands=c("lhs", "tdopt"),
           method=c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN", "optimize"),
           ...)
  {
    ## lhs should verify that the rect makes sense
    rect <- matrix(rect, ncol=2)
    
    ## XX a predictive grid
    XX <- tgp.cands(prev$obj, NN, cands, rect, X)

    ## add optim results in as a predictive location
    XX <- rbind(XX, as.numeric(prev$progress[1,1:nrow(rect)]))
    Xboth <- rbind(X,XX)
      
    ## fit a tgp model
    out <- model(X=X, Z=Z, XX=XX, improv=improv, ...)
      
    ## find the predicted minimum
    m <- which.min(c(out$Zp.mean, out$ZZ.mean))
    Xm <- Xboth[m,]
      
    ## find the optimum with kriging, and record in opt
    opt <- optim.ptgpf(Xm, rect, out, method)
    opt <- data.frame(matrix(c(opt$par, opt$value), nrow=1))
    names(opt) <- c(paste("x", 1:nrow(rect), sep=""), "z")
    
    ## X & from tgp-improv
    ir <- out$improv[,2]
    Ximprov <- matrix(XX[ir <= improv[2],], nrow=sum(ir <= improv[2]))

    ## assemble the info about the current minimum, and return
    as <- data.frame(improv=max(out$improv[,1]))
    r <- list(X=Ximprov, progress=cbind(opt, as), obj=out)
    return(r)
  }
