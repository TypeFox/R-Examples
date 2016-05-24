# The EMbC Package for R
# 
# Copyright 2013, 2014, 2015 Joan Garriga <jgarriga@ceab.csic.es>, Aitana Oltra <aoltra@ceab.csic.es>, John R.B. Palmer <johnrbpalmer@gmail.com>, Frederic Bartumeus <fbartu@ceab.csic.es>
#   
# EMbC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
# 
# EMbC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses.


# Constructors
# ------------

#' General pourpose multivariate binary Clustering (EMbC)
#' 
#' \code{embc} implements the core function of the Expectation-Maximization multivariate binary clustering.
#' 
#' @param X The input data set. A multivariate matrix where each row is a data point and each column is an input feature (a variable).
#'
#' @param U A multivariate matrix with same dimension as X with the values of certainty associated to each corresponding value in X. Ceartainties assign reliability to the data points so that the less reliable is a data point the less its leverage in the clustering. By default certainties are set to one (no uncertainty in any value in X). 
#'
#' @param stdv a vector with bounds for the maximum precission of clusters, given as minimum standard deviation for each variable, (by default is set to rep(1e-30,ncol(X))
#'
#' @param maxItr A limit to the number of iterations in case of slow convergence (defaults to 200).
#'
#' @param info Level of information shown at each step:
#' info=0 (default) shows step likelihood, number of clusters, and number of changing labels;
#' info=1 includes clustering statistics;
#' info=2 includes delimiters information;
#' info<0 supresses any step information.
#'
#' @return Returns a binClst object.
#'  
#' @export
#' @rdname embc
#'
#' @examples
#'
#' # -- apply EMbC to the example set of data points x2d ---
#' mybc <- embc(x2d@@D)

embc <- function(X,U=NULL,stdv=NULL,maxItr=200,info=0){
	if (is.null(U)) U <- matrix(rep(1,length(X)),dim(X))
	if (is.null(stdv)) stdv <- rep(1e-30,dim(X)[2])
	bC <- new('binClst',X=X,U=U,stdv=stdv)
  colnames(bC@X) <- paste('X',seq(dim(bC@X)[2]),sep='')
	colnames(bC@U) <- paste('U',seq(dim(bC@X)[2]),sep='')
	return(clst(bC,maxItr,info))
	}

# validity function for stbc constructors
checkConstructorCall <- function(obj,stdv,scv){
	errors <- character()
	if (nrow(obj)==0){
		msg <- "path length is ZERO !?"
		errors <- c(errors, msg)
		}
	if (any(is.na(obj))){
		msg <- "NA values in the path !?"
		errors <- c(errors, msg)
		}
	if (length(stdv)>0 & length(stdv)!=2){
		msg <- "invalid stdv specification"
		errors <- c(errors, msg)
		}
	if (!(scv %in% c('None','height','azimuth','rheight','rheight2','rheight3'))){
		msg <- "solar covariate must be either 'height' or 'azimuth'"
		errors <- c(errors, msg)
		}
	if (length(errors) != 0) cat(errors,'\n')
	return(errors)
	}

#' speed/turn bivariate binary Clustering.
#' 
#' \code{stbc} is a specific constructor for movement ecology pourposes. By default it implements a bivariate (speed/turn) clustering for behavioural annotation of animals' movement trajectories. Alternatively, it can perform a trivariate clustering by including the solar position covariate (i.e. solar height or solar azimuth) as a daytime indicator.
#'
#' @param obj
#'
#' A \code{data.frame} object with (timeStamp,lon,lat) values in columns 1:3 respectively. Timestamps must be given
#' as.POSIXct() with specific format "\%Y-\%m-\%d \%H:\%M:\%S". Further columns of associated data are allowed and will be included in the \link{binClstPath_instance} @@pth slot.
#'
#' A \code{Move} object from the "move" R-package.
#'
#' A \code{list} of trajectories given either as \code{data.frame} or \code{Move} objects, to perform a joined clustering of all of them. This is mainly intended to perform analysis at population level.
#'
#' @param spdLim A speed limit for automatic detection of outliers. Trajectory locations with associated values of speed above
#' the spdLim are not eliminated but will play no part in the clustering. By default is set to 40 m/s.
#'
#' @param smth A smoothing time interval in hours. This is used to estimate local values of speed and turn computed as an average over a time window centered at each location.
#'
#' @param scv A solar position covariate to be used as a daytime indicator. It can be either 'height' (the solar height in degrees above the horizon) or 'azimuth' (the solar azimuth in degrees from north). If it is used, a trivarate clustering is performed, increasing to a maximum of 8 the number of clusters (behaviours) that can potentially be identified. By default this value is set to None (i.e. perform the standard bivariate speed/turn clustering). 
#'
#' @param stdv a vector with bounds for the maximum precission of clusters, given as minimum standard deviation for each variable, (by default is set to 0.1 m/s for velocities and 5 degrees for turns).
#'
#' @param maxItr A limit to the number of iterations in case of slow convergence (defaults to 200).
#'
#' @param info Level of information shown at each step:
#' info=0 (default) shows step likelihood, number of clusters, and number of changing labels;
#' info=1 includes clustering statistics;
#' info=2 includes delimiters information;
#' info<0 supresses any step information.
#'
#' @return Returns a binClstPath object.
#'  
#' @export
#' @rdname stbc
#'
#' @examples
#' # -- apply EMbC to the example path --
#' mybcp <- stbc(expth)
#' # --- binary clustering of a Move object ---
#' require(move)
#' mybcm <- stbc(move(system.file("extdata","leroy.csv.gz",package="move")))
#'\dontrun{
#' # --- binary clustering of a stack of trajetories ---
#' mybcm <- stbc(list(mypth1,mypth2,mypth3))
#'}
setGeneric("stbc",
	function(obj,stdv=c(0.1,5*pi/180),spdLim=40,smth=0,scv='None',maxItr=200,info=0)
	{standardGeneric("stbc")})

#' @rdname stbc
setMethod("stbc",signature(obj="data.frame"),function(obj,stdv,spdLim,smth,scv,maxItr,info){
	errors <- checkConstructorCall(obj,stdv,scv)
	if (length(errors)==0){
		bCP <- new('binClstPath',obj,stdv,spdLim,smth,scv)
		return(clst(bCP,maxItr,info))}
	})

#' @rdname stbc
setMethod("stbc",signature(obj="Move"),function(obj,stdv,spdLim,smth,scv,maxItr,info){
	bCM <- new('binClstMove',obj,stdv,spdLim,smth,scv)
	return(clst(bCM,maxItr,info))})
	
#' @rdname stbc
setMethod("stbc",signature(obj="list"),function(obj,stdv,spdLim,smth,scv,maxItr,info){
	stck <- new('binClstStck',obj,stdv,spdLim,smth,scv)
	stck@bC <- clst(stck@bC,maxItr,info)
	iniLoc <- 1
	for (i in seq(length(stck@bCS))) {
		stck@bCS[[i]]@m <- ncol(stck@bCS[[i]]@X)
		stck@bCS[[i]]@k <- 2**stck@bCS[[i]]@m
		stck@bCS[[i]]@n <- nrow(stck@bCS[[i]]@X)
		stck@bCS[[i]]@C <- getColors(stck@bCS[[i]]@k)
		stck@bCS[[i]]@R <- stck@bC@R
		stck@bCS[[i]]@A <- stck@bC@A[iniLoc:(iniLoc+(stck@bCS[[i]]@n-1))]
		stck@bCS[[i]]@W <- stck@bC@W[iniLoc:(iniLoc+(stck@bCS[[i]]@n-1)),]
		iniLoc <- iniLoc + stck@bCS[[i]]@n
		}
	return(stck)})
