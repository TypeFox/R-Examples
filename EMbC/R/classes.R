# The EMbC Package for R
#
# Copyright 2013, 2014, 2015 Joan Garriga <jgarriga@ceab.csic.es>, Aitana Oltra <aoltra@ceab.csic.es>, John R.B. Palmer <johnrbpalmer@gmail.com>, Frederic Bartumeus <fbartu@ceab.csic.es>
#
#   EMbC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
#
# EMbC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses.

# Class: binClst
# --------------

#' @title Binary Clustering Class
#'
#' @description \code{binClst} is a generic multivariate binary clustering object.
#'
#' @slot X The input data set. A multivariate matrix where each row is a data point and each column is an input feature (a variable).
#' @slot U A multivariate matrix with same dimension as X with the values of certainty associated to each corresponding value in X. Ceartainties assign reliability to the data points so that the less reliable is a data point the less its leverage in the clustering. By default certainties are set to one for all variables of all data points.
#' @slot stdv A numeric vector with variable specific values for minimum standard deviation.
#' @slot m The number of input features.
#' @slot k The number of clusters.
#' @slot n The number of observations (data points).
#' @slot R A matrix with the values delimiting each binary region (the \code{Reference} values).
#' @slot P A list with the GMM (Gaussian Mixture Model) parameters. Each element of the list corresponds to a component of the GMM and it is a named-sublist itself, with elements '$M' (the component's mean) and '$S' (the component's covariance matrix).
#' @slot W A n*k matrix with the likelihood weights.
#' @slot A A numeric vector with the clustering labels (annotations) for each data-point (the basic output data). Labels are assigned based on the likelihood weights. Only in case of equal likelihoods the delimiters are used as a further criterion to assign labels.
#' @slot L The values of likelihood at each step of the optimization process.
#' @slot C Default color palette used for the plots. Can be changed by means of the setc() function.
#'
setClass("binClst",
	representation(
		X="matrix",
		U="matrix",
		stdv="numeric",
		m="integer",
		k="numeric",
		n="integer",
		R="matrix",
		P="list",
		W="matrix",
		A="numeric",
		L="numeric",
		C="character"),
	validity=function(object) {
		errors <- character()
		if (any(dim(object@X)!=dim(object@U))){
			msg <- "X,U have different dimension"
			errors <- c(errors, msg)
		}
		if (dim(object@X)[2]!=length(object@stdv)){
			msg <- "invalid stdv specification"
			errors <- c(errors, msg)
		}
		if (length(which(object@U<0|object@U>1))>0){
			msg <- "Uncertainty values must be in range (0:1)"
			errors <- c(errors, msg)
		}
		if (any(is.na(object@X))){
			msg <- "NA/NaN values found"
			errors <- c(errors, msg)
		}
		if (length(errors) == 0) TRUE else errors
	})


# Class: binClstPath
# ------------------

#' @title Binary Clustering Path Class
#'
#' @description \code{binClstPath} is a \code{binClst} subclass for fast and easy speed/turn-clustering of movement trajectories. The input trajectory is given as a data.frame with, at least, the columns (timeStamp,longitude,latitude). This format is described in detail in the class constructor \link{stbc}. As a \code{binClst} subclass, this class inherits all slots and functionality of its parent class.
#'
#' @slot pth A data.frame with the trajectory timestamps and geolocation coordinates, plus eventual extra columns that were included in the input path data frame, (see the \link{stbc} constructor).
#' @slot spn A numeric vector with the time intervals between locations (in seconds).
#' @slot dst A numeric vector with the distances between locations (in meters). We use loxodromic computations.
#' @slot hdg A numeric vector with local heading directions (in radians from North). We use loxodromic computations.
#' @slot bursted A logical value indicating whether the \code{binClstPath} instance has already been bursted. As bursting can be computationally demanding for long trajectories, an instance is bursted only when a burst wise representation of the trajectory' is requested for the first time, (unless this value is changed to FALSE).
#' @slot tracks If bursted=TRUE, a \code{SpatialLinesDataFrame} object ("sp" R-package) with the bursted track segments.
#' @slot midPoints If bursted=TRUE, a \code{SpatialPointsDataFrame} object ("sp" R-package) with the bursted track midpoints.
#'
setClass("binClstPath",
	representation(
		pth="data.frame",
		spn="numeric",
		dst="numeric",
		hdg="numeric",
		bursted="logical",
		tracks="SpatialLinesDataFrame",
		midPoints="SpatialPointsDataFrame"),
	contains=c("binClst"))

setMethod("initialize","binClstPath",function(.Object,pth=data.frame(),stdv,spdLim,smth,scv){
	if (nrow(pth)>0){
		.Object@pth <- bCPStd(pth)
		.Object@spn <- c(spanTime(.Object@pth),0)
		.Object@dst <- c(loxDst(.Object@pth),0)
		.Object@hdg <- c(loxTht(.Object@pth),0)
		if (scv!='None'){
			.Object@X <- cbind(getSolarPos(.Object@pth,scv),getSpeed(.Object),getTurns(.Object))
      colnames(.Object@X) <- c(scv,'velocity (m/s)','turn (rad)')
			.Object@stdv <- c(c(5.0,5.0,0.05,0.01,0.005)[which(c('azimuth','height','rheight','rheight2','rheight3')==scv)],stdv)
			.Object@U <- cbind(rep(1,nrow(pth)),stCertainty(.Object))
			colnames(.Object@U) <- c('scv.Certainty','vlc.Certainty','trn.Certainty')
			if (spdLim>0){
				.Object@U[which(.Object@X[,2]>spdLim),] <- c(0,0,0)
				.Object@X[which(.Object@X[,2]>spdLim),2] <- 0}
		} else {
			.Object@X <- cbind(getSpeed(.Object),getTurns(.Object))
			colnames(.Object@X) <- c('velocity (m/s)','turn (rad)')
			.Object@stdv <- stdv
			.Object@U <- stCertainty(.Object)
			colnames(.Object@U) <- c('vlc.Certainty','trn.Certainty')
			if (spdLim>0){
				.Object@U[which(.Object@X[,1]>spdLim),] <- c(0,0)
				.Object@X[which(.Object@X[,1]>spdLim),1] <- 0}
		}
		.Object@bursted <- FALSE
		if (smth>0) .Object <- priorSmth(.Object,smth)
	}
	.Object})


# Class: binClstMove
# ------------------

#' @title Binary Clustering Move Class
#'
#' @description \code{binClstMove} is a \code{binClstPath} subclass for speed/turn-clustering of \code{Move} objects from the \code{move} R-package. This class inherits all slots and functionality of \code{binClstPath} and \code{Move} objects.
#'
setClass("binClstMove",
	contains=c("binClstPath","Move"))

setMethod("initialize","binClstMove",function(.Object,obj,stdv,spdLim,smth,scv){
	bCP <- new('binClstPath',data.frame(obj$study.local.timestamp,obj@coords),stdv,spdLim,smth,scv)
	for(i in slotNames(bCP))
		slot(.Object,i) <- slot(bCP,i)
	for(i in slotNames(obj))
		slot(.Object,i)<-slot(obj,i)
	.Object})


# Class: binClstStck
# ----------------------

#' @title Binary Clustering Stack Class
#'
#' @description \code{binClstStck} is a special class for population level speed/turn-clustering of movement trajectories, given either as path data.frames or \code{move} objects.
#'
#' @slot bCS A list of either \code{binClstPath} or \code{binClstMove} objects, depending on how the input paths are given.
#' @slot bC A \code{binClst} instance with the global speed/turn clustering of the paths in the stack.
#'
setClass("binClstStck",
	representation(
		bCS="list",
		bC="binClst")
	)

setMethod("initialize","binClstStck",function(.Object,obj,stdv,spdLim,smth,scv){
	.Object@bCS <- lapply(obj,function(pth) {
		if (class(pth)=='data.frame') new('binClstPath',pth,stdv,spdLim,smth,scv)
		else if (class(pth)=='Move') new('binClstMove',pth,stdv,spdLim,smth,scv)
		})
	X <- .Object@bCS[[1]]@X
	U <- .Object@bCS[[1]]@U
	for (i in 2:length(.Object@bCS)){
		X <- rbind(X,.Object@bCS[[i]]@X)
		U <- rbind(U,.Object@bCS[[i]]@U)
		}
	if (scv!='None')
		.Object@bC <- new('binClst',X=X,U=U,stdv=c(5.0,stdv))
	else
		.Object@bC <- new('binClst',X=X,U=U,stdv=stdv)
	.Object})
