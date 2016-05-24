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
MLESpatialProcess <- function(x, y, theta.grid=NULL, par.grid=NULL, lambda.grid=NULL, 
                             cov.function = "stationary.cov", 
                             cov.args = list(Covariance = "Matern", smoothness = 1), 
                             optim.args = NULL, ngrid = 10, niter = 15, tol = 0.01, 
                             Distance = "rdist", nstep.cv = 50, verbose = FALSE, 
                             doMKrig=FALSE, ...)
{
  if(!doMKrig) {
    # if using Krig not mKrig
    
    if(!is.null(par.grid))
      stop("par.grid argument not supported when using MLESpatialProcess with Krig instead of mKrig")
    
    # save arguments for Krig as a list
    
    KrigCallingList <- c(
      # list to pass to the objective function
      info <- list( x = x, 
                    Y = y, 
                    cov.function = cov.function,
                    cov.args = cov.args, 
                    method = "REML",
                    nstep.cv = nstep.cv,
                    give.warnings = FALSE, 
                    Distance = Distance),
      list( ...)
    )
    # objective function used for grid search and optimization	
    objective.fn <- function(ltheta, returnAll = FALSE) {	
      thetaList<- list( theta=exp( ltheta))
      hold <- do.call(Krig, 
                      c( KrigCallingList,thetaList) )[c("lambda.est", "gcv.grid")]
      minus.lPLike <- hold$lambda.est["REML", 5]
      # add this evalution to an  object
      # (i.e. here a matrix) in the calling frame
      withAddedRow <- rbind(
        get("capture.evaluations", envir = capture.env),
        c(exp(ltheta),  hold$lambda.est["REML",  ])
      )
      if( verbose){
        print( c(exp(ltheta),  hold$lambda.est["REML",  ]))
      }
      assign("capture.evaluations", withAddedRow, envir = capture.env)
      # return all serach information or just profile like value.        
      if (returnAll) {
        return(hold)
      } else {
        return(minus.lPLike)
      }
    }	
    # set up matrix to cature all  evaluations from within optimization
    capture.evaluations <- matrix(NA, ncol = 7, nrow = 1, 
                                  dimnames = list(NULL, 
                                                  c("theta", "lambda.MLE", 
                                                    "trA", "GCV", "sigmaGCV", 
                                                    "lnProfileLike","warning")))
    capture.env <- environment()
    #
    # if grid for ranges is missing use  quantiles of pairwise
    #distances among data.
    #
    if (is.null(theta.grid)) {
      # Distances between locations
      pairwiseD<- get(Distance)(x,x)
      pairwiseD<- pairwiseD[col(pairwiseD) > row( pairwiseD) ]
      theta.range <- quantile(pairwiseD , c(0.03, 0.97))
      theta.grid <- seq(theta.range[1], theta.range[2], , ngrid)
    }
    if (length(theta.grid) == 2) {
      theta.grid <- seq(theta.grid[1], theta.grid[2], , ngrid)
    }
    ngrid <- length(theta.grid)
    minus.REML <- rep(NA, ngrid)
    # object to hold log likelihood surface as a function of log lambda and 	theta
    logLikelihoodSurface<- list( x = matrix( theta.grid, nrow= ngrid, ncol=nstep.cv ),
                                 y = matrix(NA,nrow=ngrid, ncol= nstep.cv ),
                                 z = matrix( NA, nrow= ngrid,  ncol=nstep.cv))	                       
    # grid search over theta likelihood maximized over
    # lambda (sill and nugget) for each theta
    for (j in 1:ngrid) {
      out<- objective.fn(log(theta.grid[j]),  returnAll=TRUE)
      minus.REML[j] <- out$lambda.est["REML", 5]
      logLikelihoodSurface$y[j,]<- log( out$gcv.grid[,1])
      logLikelihoodSurface$z[j,]<- -1*out$gcv.grid[,7]
    }
    #	temp <- cbind(theta.grid, -minus.REML)
    #	dimnames(temp) <- list(NULL, c("theta", "logProfileLike"))
    # best point for theta from grid search
    # NOTE: due to past convenion in Krig  - log likelihood is computed
    # so this quantity is _minimized_
    IMIN <- which.min( minus.REML)
    if (IMIN == 1 | IMIN == ngrid) {
      warning("theta (range parameter) at end of search interval:", fill = TRUE)
      theta.MLE <- theta.grid[IMIN]
      REML.MLE<- -1* minus.REML[IMIN]
    } else {
      # starting interval  for  1-d  optimization
      lstart <- log(theta.grid)[c(IMIN-1, IMIN+1) ]
      out <- optimize(f = objective.fn, interval = c(lstart[1], lstart[2]), maximum = FALSE)
      theta.MLE <- exp(out$minimum)
      REML.MLE <- -1 * out$objective
    }
    eval.grid <- capture.evaluations[-1, ]
    eval.grid[, 6] <- -1 * eval.grid[, 6]
    # sort on theta!
    ind<- order( eval.grid[,1])
    eval.grid<- eval.grid[ind,]
    #		
    thetaList<- list( theta=theta.MLE)
    hold <- do.call(Krig, c( KrigCallingList, thetaList) )[c("lambda.est", "rho.MLE", "shat.MLE")]
    pars<- c( theta.MLE, hold$lambda.est[6,1], hold$rho.MLE, hold$shat.MLE)
    names( pars) <- c( "theta", "lambda", "rho", "sigma")
    return(list(
      pars = pars,
      logLikelihood = REML.MLE,
      eval.grid = eval.grid,
      logLikelihoodSurface = logLikelihoodSurface,
      lambda.est = hold$lambda.est,
      call = match.call() )
    )
  }
  else{
    #if doMKrig
    
    #theta.grid should only be used with Krig, not mKrig
    if(!is.null(theta.grid))
      stop("theta.grid argument not supported when using MLESpatialProcess with mKrig")
    
    #set lambda.grid if unset to be on logarithmic scale from 10^-6 to 10
    if(is.null(lambda.grid))
      lambda.grid = 10^seq(-6, 1, by=1)
    else if(length(lambda.grid) == 2)
      lambda.grid = 2^seq(log(lambda.grid[1], 2), log(lambda.grid[2], 2), length=ngrid)
    
    # check if any parameters in par.grid have only two elements.  If so, 
    # expand the grid
    if(!is.null(par.grid)) {
      for(p in 1:length(par.grid)) {
        vals = par.grid[[p]]
        if(length(vals) == 2)
          par.grid[[p]] = seq(log(vals[1], 2), log(vals[2], 2), length=ngrid)
      }
    }
    par.grid <- make.surface.grid(par.grid)
    
    #if too many parameters to try, give warning
    if(nrow(par.grid) == 0)
      par.grid.length = 1
    else
      par.grid.length = nrow(par.grid) ^ ncol(par.grid)
    if(par.grid.length*length(lambda.grid) > 100)
      warning("The parameter and lambda grids are large, so the likelihood maximization process may take a long time")
    
    #list for calling mKrig.MLE with everything but lambda
    #
    if(supportsArg(cov.function, "Distance"))
      cov.args$Distance=Distance
    mKrig.callingList = c(list(x=x, y=y, cov.fun=cov.function, cov.args=cov.args), 
                          list(...))
    mKrig.MLE.callingList = c(list(par.grid=par.grid, relative.tolerance=tol, 
                                   verbose=verbose), 
                              mKrig.callingList)
    
    # ensure find.trA and lambda.profile are FALSE for MLE 
    # since there's no reason for them here
    mKrig.MLE.callingList$find.trA = FALSE
    mKrig.MLE.callingList$lambda.profile = FALSE
    
    #call mKrig.MLE for each lambda in lambda.grid
    evals = data.frame()
    full.par.grid = data.frame()
    lnLik = data.frame()
    for(l in 1:length(lambda.grid)) {
      out = do.call("mKrig.MLE", c(mKrig.MLE.callingList, lambda=lambda.grid[l]))
      evals = rbind(evals, out$summary)
      full.par.grid = rbind(full.par.grid, out$par.grid)
    }
    
    #get log-likelihood surface data
    if(nrow(full.par.grid) != 0)
      params = cbind(as.matrix(full.par.grid), evals[,6])
    else
      params = cbind(evals[,6])
    lnLik = evals[,2]
    
    #create grid.list of all parameter values to interpolate to
    paramsInterp = list()
    for(col in 1:ncol(params)) {
      colRange = range(params[,col])
      paramsInterp = c(paramsInterp, list(seq(colRange[1], colRange[2], length=60)))
    }
    paramsInterp = make.surface.grid(paramsInterp)
    
    #interpolate surface with thin-plate splines
    tps = Tps(params, lnLik)
    lnLikInterp = predict(tps, paramsInterp)
    
    #get maximum of thin-plate spline lnLik surface
    maxI = which.max(lnLikInterp)
    maxParams = paramsInterp[maxI,]
    
    #exponentiate log lambda
    maxParams[length(maxParams)] = exp(maxParams[length(maxParams)])
    
    #convert max parameter vector to a list that can be used to call covariance
    covParamsNames = names(full.par.grid)
    if(length(maxParams) == 1)
      cov.arg.list=NULL
    else
      cov.arg.list = as.list(maxParams[1:(length(maxParams)-1)])
    names(cov.arg.list) = covParamsNames
    
    #make arguments for mKrig.MLE.joint
    #
    if(is.null(optim.args))
      optim.args = list(method = "BFGS", 
                        control=list(fnscale = -1, 
                                     ndeps = rep(log(1.1), length(maxParams)), 
                                     reltol=1e-02, maxit=3))
    
    if(!is.null(tol))
      optim.args$control$reltol = tol
    
    lambda.guess = maxParams[length(maxParams)]
    cov.params.guess = cov.arg.list
    
    # perform final joint optimization with initial guess as the Tps maximum,
    # and return final jointly optimized mKrig object
     out<- do.call("mKrig.MLE.joint", c(mKrig.callingList, 
                                 list(lambda.guess=lambda.guess, 
                                      cov.params.guess=cov.params.guess, 
                                      optim.args=optim.args, 
                                      verbose=verbose, 
                                      find.trA.MLE=FALSE)))
    # replace strange call that is generated form this internal 
    # call to mKrig with the top level one
      out$call<- match.call()                               
      return( out)                                
    
#     # perform mKrig at Tps maximum parameters
#     mKrig.callingList$cov.args = c(cov.args, cov.arg.list)
#     mKrig.callingList$lambda = maxParams[length(maxParams)]
#     return(do.call("mKrig", mKrig.callingList))
  }
}




