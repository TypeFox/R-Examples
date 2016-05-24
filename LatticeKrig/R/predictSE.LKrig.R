# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
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

"predictSE.LKrig" <- function(object, xnew = NULL, Znew = object$Z, 
	verbose = FALSE, ...) {
	if (is.null(object$Mc)) {
		stop("need to include the sparse cholesky decompostion in LKrig object \r\nin calling LKrig set return.cholesky = TRUE")
	}
	if (is.null(object$LKinfo$fixedFunction)) {
		stop("Standard errors not supported without a fixed component")
	}
	# set some local variables
	NG <- nrow(xnew)
	lambda <- object$lambda
	rho <- object$rho.MLE
	sigma2 <- lambda * rho
	weights <- object$weights
# NOTE throughout the "w" added to variables, e.g. wX, wk0, is the sqrt(weights)
# not the weights	
	LKinfo <- object$LKinfo
	if( is.null( xnew)){
		xnew<- object$x
	}
# figure out if extra covariates should be included
if (is.null(Znew) & (object$nZ > 0)) {
		Znew <- object$Z
	}
	distance.type <- LKinfo$distance.type
	if(is.null(object$wU) | is.null( object$wX) ){
		stop("LKrig object must have the matrices wU and wX 
		 (see return.wXandwU = TRUE in LKrig)")}
# fixed effects matrix at the new locations		 
	t0 <- t(do.call(object$LKinfo$fixedFunction, c(list(x = xnew, 
		Z = Znew, distance.type = distance.type), LKinfo$fixedFunctionArgs)))
	Omega <- object$Omega
	#
#	wS  <- diag.spam(sqrt(weights))
#	wk0 <- wS%*% LKrig.cov(object$x, xnew, LKinfo = LKinfo)
	wk0 <- LKrigCovWeightedObs( xnew, object$wX, LKinfo = LKinfo)
    # use the predict computations but with covariance vector
    # as "observations"
    # it may not be obvious why this formula makes sense!	
	hold <- LKrig.coef(
	                   Mc = object$Mc,
	                   wX = object$wX, 
		               wU = object$wU,
		               wy = wk0,
		           lambda = lambda,
		          verbose = verbose)
	# version of 'c'coefficents for usual Kriging model
	c.mKrig <- (wk0 - object$wU %*% hold$d.coef
	                        - object$wX %*% hold$c.coef)/lambda
	d.coef <- hold$d.coef
	# colSums used to find multiple quadratic forms  
	#e.g.  diag(t(x) %*%A%*%x) == colSums( x* (A%*%(x)))
	temp1 <- rho * (colSums(t0 * (Omega %*% t0)) - 2 * colSums(t0 * 
		d.coef))
	# find marginal variances -- trival in the stationary case!
	temp0 <- rho * (LKrig.cov(xnew, LKinfo = LKinfo, marginal = TRUE) - 
		colSums(wk0 * c.mKrig))
	temp <- temp0 + temp1
	return(sqrt(temp))
}
