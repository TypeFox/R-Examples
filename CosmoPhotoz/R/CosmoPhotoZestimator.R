#  R package CosmoPhotoz file R/CosmoPhotoZestimator.R
#  Copyright (C) 2014 COIN
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

#' @title Photometric redshift  estimation from a training dataset and a test dataset
#'
#' @description \code{CosmoPhotoZestimator} returns photometric redshift estimated from photometric 
#' data and a training dataset with photometry and spectroscopy. The estimation is 
#' based on generalized linear models (see \code{\link{glmTrainPhotoZ}} and \code{\link{glmPredictPhotoZ}}).
#' 
#' @import shiny
#' @param trainData vector containing spectroscopic redshift data and photometry (at least one column shall be called redshift)
#' @param testData vector containing spectroscopic redshift data and photometry (at least one column shall be called redshift)
#' @param numberOfPcs an integer indicating the number of principal components to consider
#' @param method a string containing the chosen GLM method. Two options are available: \code{Frequentist} will use the function  \code{\link{glm}} from the package \code{stats}; \code{Bayesian} will use the function \code{\link{bayesglm}} from the package \code{arm}
#' @param family a string containing \code{gamma} or \code{inverse.gaussian} (a description of the error distribution and link function to be used in the model)
#' @param robust a boolean indicating if robust PCA should be used or not
#' @return a vector with the estimated photometric redshifts
#' @examples
#' \dontrun{
#' # Load the data
#' data(PHAT0train)
#' data(PHAT0test)
#' 
#' # Run the analysis
#' photoZest <- CosmoPhotoZestimator(PHAT0train, PHAT0test, 6)
#' 
#' # Create a boxplot showing the results
#' plotDiagPhotoZ(photoz = photoZest, specz = PHAT0test$redshift, type = "box")
#' }
#' 
#' @usage CosmoPhotoZestimator(trainData, testData, numberOfPcs, method, family, robust)
#' 
#' @author Alberto Krone-Martins, Rafael S. de Souza
#' 
#' @keywords utilities
#' @export
CosmoPhotoZestimator <- function(trainData, testData, numberOfPcs=4, method="Bayesian", family="gamma", robust=TRUE) {
  # Combine the training and test data and calculate the principal components
  redshift <- NULL # <- this is just to prevent a NOTE from CRAN checks
  PC_comb <- computeCombPCA(subset(trainData, select=c(-redshift)),
                            subset(testData,  select=c(-redshift)),
                            robust=robust)
  Trainpc <- cbind(PC_comb$x, redshift=trainData$redshift)
  Testpc <- PC_comb$y

  # Dynamic generation of the formula based on the user selected number of PCs
  formMa <- "poly(Comp.1,2)*poly(Comp.2,2)*"
  formMb <- paste(names(PC_comb$x[3:numberOfPcs]), collapse="*")
  formM <- paste("redshift~",formMa, formMb, sep="")

  # Fitting
  Fit <- glmTrainPhotoZ(Trainpc, formula=eval(parse(text=formM)), 
                        method=method, family=family)

  # Photo-z estimation
  photoz <- predict(Fit$glmfit, newdata=Testpc, type="response")
  
  return(photoz)
}
