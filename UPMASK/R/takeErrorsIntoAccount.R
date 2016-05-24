#  R package UPMASK file R/takeErrorsIntoAccount.R
#  Copyright (C) 2014 Alberto Krone-Martins
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

#' @title Take Errors Into Account for UPMASK analysis
#' 
#' @description Based on a data frame containing measurements and errors, the 
#' \code{takeErrorsIntoAccount} will produce another data frame where each measurement of 
#' the original data frame is replaced by another value taken from a random distribution. 
#' The implemented error model is gaussian, so each value of the output data frame will
#' be a random sampling from a gaussian distribution where the mean is the value in the 
#' original data frame (indicated by the \code{photometricDataIndexes} column argument) 
#' and the standard deviation is the value from its corresponding error (indicated by the 
#' \code{photometricErrorDataIndexes} column argument). The newly constructed dataframe
#' is returned by the function.
#' 
#' The user can adapt this function so it can take any error model into account during the
#' UPMASK analysis.
#' 
#' @param originalData a data frame to use as the baseline
#' @param dataIndexes an array of integers indicating the columns corresponding to the measurements
#' @param errorIndexes an array of integers indicating the columns corresponding to the errors
#' 
#' @return A data frame with the new values sampled from the error distributions.
#'
#' @references \href{http://dx.doi.org/10.1051/0004-6361/201321143}{Krone-Martins, A. & Moitinho, A., A&A, v.561, p.A57, 2014}
#'
#' @examples
#' # Create a simple data set with the values and errors
#' toyDataDF <- data.frame(x=runif(10, 0, 10), dx=rep(0.2, 10), y=runif(10, 0, 10), 
#'                         dy=rep(0.1, 10))
#' 
#' # Apply the error models to create another data frame
#' newToyDataDF <- takeErrorsIntoAccount(toyDataDF, c(1,3), c(2,4))
#' 
#' # Plot the results
#' plot(toyDataDF$x, toyDataDF$y)
#' points(newToyDataDF$x, newToyDataDF$y, pch=19, cex=0.8, col="red")
#' 
#' # Clean the environment
#' rm(list=c("toyDataDF", "newToyDataDF"))
#'  
#' @usage takeErrorsIntoAccount(originalData, dataIndexes, errorIndexes)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords utilities
#' @export
#
takeErrorsIntoAccount <- function(originalData, dataIndexes, errorIndexes) {
  
  outDF <- originalData
  
  # The default is a simple gaussian error model
  # The user can modify this function as needed
  for(i in 1:length(dataIndexes)) {
  	outDF[,dataIndexes[i]] <- rnorm(length(outDF[,1]), 
  	                                             outDF[,dataIndexes[i]], 
  	                                             outDF[,errorIndexes[i]])
  }
  	
  # That is all folks!
  return(outDF)
}