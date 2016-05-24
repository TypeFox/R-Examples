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
MLESpatialProcess.fast <- function(x, y, lambda.start=NULL,
                                   theta.start = NULL, 
     cov.function = "stationary.cov", 
	     cov.args = list(Covariance = "Matern", smoothness = 1), 
	 relative.tolerance = 1e-3, Distance = "rdist", 
	verbose=FALSE,
	...) {
# save arguments for Krig as a list
	
	mKrigCallingList <- c(
	# list to pass to the objective function
    info <- list( x = x, 
                  y = y, 
       cov.function = cov.function,
           cov.args = cov.args, 
     	   Distance = Distance),
     	    list( ...)
    	   )
# objective function used for grid search and optimization	
	objective.fn <- function(par, returnTrA=FALSE) {	
            parList<- list( theta=exp(par[1]), lambda = exp(par[2]),
                           find.trA = returnTrA )
		hold <- do.call(mKrig, 
		     c( mKrigCallingList,parList) )[c("lambda.fixed", 
                                                      "rho.MLE.FULL",
                                                      "eff.df", "GCV",
                                                      "sigma.MLE.FULL",
                                                      "lnProfileLike.FULL")
                                                    ]
            hold<- unlist( hold)
		logPLike <- hold[6]
		# add this evalution to an  object
		# (i.e. here a matrix) in the calling frame
	
        withAddedRow <- rbind(
             get("capture.evaluations", envir = capture.env),
             c( par[1], hold )
             )
             if( verbose){
             	print( c(par, hold ) )
             }
		assign("capture.evaluations", withAddedRow, envir = capture.env)
		# return all search information or just profile like value.        
			return(logPLike)
	}	
	# set up matrix to cature all  evaluations from within optimization
	capture.evaluations <- matrix(NA, ncol = 7, nrow = 1, dimnames = list(NULL, 
		c("theta", "lambda", "rhoMLE", "trA", "GCV", "sigmaMLE", "lnProfileLike")))
	capture.env <- environment()
	#
	# if grid for ranges is missing use  quantiles of pairwise
#distances among data.
#
if (is.null(theta.start)) {
		# Distances between locations
		pairwiseD<- get(Distance)(x,x)
		pairwiseD<- pairwiseD[col(pairwiseD) > row( pairwiseD) ]
		theta.range <- quantile(pairwiseD , c(0.03, 0.97))
		theta.start <- median( pairwiseD)
	}
if( is.null( lambda.start)){
	lambda.start<- .5	
}	
#
# NOTE: due to past convenion in Krig  - log likelihood is computed
# so this quantity is _minimized_
     look <- optim(c(log(lambda.start), log(theta.start)), objective.fn, method = "BFGS", 
                control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                  ndeps = c(0.05,0.05), reltol = relative.tolerance))
#return(look)
		theta.MLE <- exp(look$par[1])
		lambda.MLE <- exp( look$par[2])
	eval.grid <- capture.evaluations[-1, ]
	eval.grid[,1]<- exp(eval.grid[,1])
	 logLike<- look$value
# sort on theta!
	ind<- order( eval.grid[,1])
	eval.grid<- eval.grid[ind,]
#
ind<- theta.MLE == eval.grid[,1]
pars <- 	 eval.grid[ind,c(1:3, 6)]
names( pars) <- c("theta", "lambda", "rho", "sigma")
		return(list(
		          pars =  pars,
		 logLikelihood = logLike,
		     eval.grid = eval.grid,
             converge  = c(look$convergence, look$counts),
		          call = match.call() )
		       )   		     
}
