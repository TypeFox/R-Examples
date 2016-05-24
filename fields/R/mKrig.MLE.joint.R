# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
mKrig.MLE.joint <- function(x, y, weights = rep(1, nrow(x)), 
                            lambda.guess = 1, cov.params.guess=NULL, 
                            cov.fun="stationary.cov", cov.args=NULL, 
                            Z = NULL, optim.args=NULL, find.trA.MLE = FALSE, 
                            ..., verbose = FALSE) {
  
  #set default optim.args if necessary
  if(is.null(optim.args))
    optim.args = list(method = "BFGS", 
             control=list(fnscale = -1, 
                          ndeps = rep(log(1.1), length(cov.params.guess)+1), 
                          reltol=1e-04, maxit=10))
  
  #check which optimization options the covariance function supports
  supportsDistMat = supportsArg(cov.fun, "distMat")
  
  #precompute distance matrix if possible so it only needs to be computed once
  if(supportsDistMat) {
    
    #Get distance function and arguments if available
    #
    Dist.fun= c(cov.args, list(...))$Distance
    Dist.args=c(cov.args, list(...))$Dist.args
    
    #If user left all distance settings NULL, use rdist with compact option.
    #Use rdist function by default in general.
    #
    if(is.null(Dist.fun)) {
      Dist.fun = "rdist"
      if(is.null(Dist.args))
        Dist.args = list(compact=TRUE)
    }
    
    distMat = do.call(Dist.fun, c(list(x), Dist.args))
  }
  
  #set cov.args for optimal performance if possible
  if(supportsDistMat)
    cov.args = c(cov.args, list(distMat=distMat, onlyUpper=TRUE))
  
  # these are all the arguments needed to call mKrig except lambda and cov.args
  mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z, cov.fun=cov.fun), 
                  list(...))
  mKrig.args$find.trA = find.trA.MLE
  
  # output matrix to summarize results
  ncolSummary = 8 + length(cov.params.guess)
  summary <- matrix(NA, nrow = 1, ncol = ncolSummary)
  dimnames(summary) <- list(NULL, c("EffDf", "lnProfLike", "GCV", "sigma.MLE", 
                                    "rho.MLE", "llambda.MLE", names(cov.params.guess), 
                                    "counts eval","counts grad"))
  
  # Define the objective function as a tricksy call to mKrig
  # if Y is a matrix of replicated data sets use the log likelihood for the complete data sets
  lnProfileLike.max <- -Inf
  temp.fn <- function(parameters) {
    # Seperate lambda from covariance parameters.
    # Optimization is over log-scale so exponentiate log-parameters.
    lambda = exp(parameters[1])
    if(length(parameters) > 1) {
      otherParams = as.list(exp(parameters[2:length(parameters)]))
      names(otherParams) = names(cov.params.guess)
    }
    else
      otherParams = NULL
    
    #get all this eval's covariance arguments using the input parameters
    cov.args.temp = c(cov.args, otherParams)
    
    # NOTE: FULL refers to estimates collapsed across the replicates if Y is a matrix
    # assign to hold the last mKrig object
    hold <- do.call("mKrig", c(mKrig.args, list(lambda = lambda),
                               cov.args.temp))
    
    #save best mKrig object to global environment
    if(hold$lnProfileLike.FULL > lnProfileLike.max) {
      out <<- hold
      lnProfileLike.max = hold$lnProfileLike.FULL
    }
    hold = hold[c("rho.MLE.FULL", "sigma.MLE.FULL", "lnProfileLike.FULL")]
    
    # add this evalution to an object (i.e. here a matrix) in the calling frame
    temp.eval <- get("capture.evaluations")
    assign("capture.evaluations", rbind(temp.eval, c(parameters, unlist(hold))), envir = capture.env)
    return(hold$lnProfileLike.FULL)
  }
  
  #
  # optimize over covariance parameters and lambda
  
  # list of covariance arguments from par.grid with right names (some R arcania!)
  # note that this only works because 1) temp.fn will search in this frame for this object
  # par.grid has been coerced to a data frame so one has a concept of a row subscript.
  
  # set up matrix to store evaluations from within optim
  capture.evaluations <- matrix(NA, ncol = 4+length(cov.params.guess), nrow = 1,
                                dimnames = list(NULL, c("lambda", names(cov.params.guess), "rho.MLE",
                                                        "sigma.MLE", "lnProfileLike.FULL")))
  capture.env <- environment()
  
  # call to optim with initial guess (on log-scale)
  init.guess = log(unlist(c(lambda.guess, cov.params.guess)))
  look <- do.call(optim, c(list(par=init.guess), list(temp.fn), optim.args))
  
  #get optim results
  optim.counts <- look$counts
  llambda.opt <- look$par[1]
  lambda.opt <- exp(llambda.opt)
  if(length(look$par) > 1) {
    params.opt <- exp(look$par[2:length(look$par)])
    params.opt <- as.list(params.opt)
    names(params.opt) <- names(cov.params.guess)
  }
  else
    params.opt=NULL
  
  # call to 1-d search
  #            opt.summary     <- optimize(temp.fn, interval= llambda.start + c(-8,8), maximum=TRUE)
  #            llambda.opt <- opt.summary$maximum
  #            optim.counts<- c(nrow(capture.evaluations)-1, NA)
  # accumulate the new matrix of lnlambda and ln likelihoods (omitting first row of NAs)
  lnLik.eval <- capture.evaluations[-1,]
  
  #exponentiate lambda and covariance parameters in lnLik.eval
  lnLik.eval[, 1:length(look$par)] = exp(lnLik.eval[, 1:length(look$par)])
  
  # calculate trace of best mKrig object if necessary
  #
  find.trA = list(...)$find.trA
  if(is.null(find.trA) || find.trA) {
    
    #get arguments for mKrig.trace
    iseed = list(...)$iseed
    NtrA = list(...)$NtrA
    
    #set iseed and NtrA to default values of mKrig if NULL
    if(is.null(iseed))
      iseed = 123
    if(is.null(NtrA))
      NtrA = 20
    
    #update best mKrig object with trace results
    out2 <- mKrig.trace(out, iseed, NtrA)
    out$eff.df <- out2$eff.df
    out$trA.info <- out2$trA.info
    np <- nrow(x)
    out$GCV <- (sum(out$residuals^2)/np)/(1 - out2$eff.df/np)^2
    if (NtrA < np)
      out$GCV.info <- (sum(out$residuals^2)/np)/(1 - out2$trA.info/np)^2
    else
      out$GCV.info <- NA
  }
  
  # save results of the best covariance model evaluation in a neat table
  summary[1, 1:ncolSummary] <- unlist(c(out$eff.df, out$lnProfileLike.FULL, 
                                 out$GCV, out$sigma.MLE.FULL, out$rho.MLE.FULL, llambda.opt, 
                                 params.opt, optim.counts))
  if (verbose) {
    cat("Summary: ", 1, summary[1, 1:ncolSummary], fill = TRUE)
  }
  
  #add summary table to output mKrig object and ensure it is still 
  #of class mKrig
  out = c(out, list(summary=summary, lnLik.eval=lnLik.eval))
  class(out) = "mKrig"
  
  return(out)
}
