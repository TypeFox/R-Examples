#  R package CosmoPhotoz file R/glmPredictPhotoZ.R
#  Copyright (C) 2014  COIN
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#' @title Predict photometric redshifts using a given glm fit object
#'
#' @description \code{glmPredictPhotoZ} computes a list of simple summary 
#' statistics for the photometric redshift estimation.
#' 
#' @param data a data.frame containing the data one wished to compute the redshift
#' @param train a trained glm object containing the fit of the model
#' @return list containing the results of the redshift estimation
#' @examples
#' \dontrun{
#' # Load the data
#' data(PHAT0train)
#' data(PHAT0test)
#' 
#' # Combine the training and test data and calculate the principal components
#' PC_comb <- computeCombPCA(subset(PHAT0train, select=c(-redshift)), 
#'            subset(PHAT0test, select=c(-redshift)),
#'            robust=FALSE) # robust is false here just to make it faster
#' Trainpc <- cbind(PC_comb$x, redshift=PHAT0train$redshift)
#' Testpc <- PC_comb$y
#' 
#' # Fitting
#' Fit <- glmTrainPhotoZ(Trainpc, formula=redshift~poly(Comp.1,2)*
#'            poly(Comp.2,2)*Comp.3*Comp.4*Comp.5*Comp.6, 
#'            method="Bayesian", family="gamma")
#' 
#' # Perform the photo-z estimation using the glmPredictPhotoZ function
#' photoz <- glmPredictPhotoZ(data=Testpc, train=Fit$glmfit)
#' specz <- PHAT0test$redshift
#' 
#' # Show a plot with the results
#' plotDiagPhotoZ(photoz$photoz, specz, "box")
#' }
#
#' @usage glmPredictPhotoZ(data, train)
#' 
#' @author Rafael S. de Souza, Alberto Krone-Martins
#' 
#' @keywords utilities
#' @export
glmPredictPhotoZ <- function(data, train){
	
 	# First some basic error control
 	if( ! is.data.frame(data) ) {
 		stop("Error in glmPredictPhotoZ :: data is not a data frame, and the code expects a data frame.")
 	}
 	
	# Now for the real work
	#Photoz<-predict(train,newdata=subset(data,select=-c(redshift)),type="response",se.fit = TRUE)
	photozObj <- predict(train, newdata=data, type="response", se.fit = TRUE)
	photoz <- photozObj$fit
	err_photoz <- photozObj$se.fit

	
	return(list(photoz=photoz, err_photoz=err_photoz))
}


