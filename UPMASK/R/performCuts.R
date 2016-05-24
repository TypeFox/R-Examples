#  R package UPMASK file R/performCuts.R
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

#' @title Perform cuts in the data
#' 
#' @description \code{performCuts} will perform cuts in the data. This function is provided
#' as a place holder, and it is empty, but it is called by \code{UPMASK}, so if the user 
#' needs to perform cuts in the data for the UPMASK analysis, this function should be 
#' tailored.
#' 
#' @param originalData a data frame to use as the baseline
#' 
#' @return A data frame.
#'
#' @examples
#' # Create a simple data set with the values and errors
#' toyDataDF <- data.frame(x=runif(10, 0, 10), dx=rep(0.2, 10), y=runif(10, 0, 10), 
#'                         dy=rep(0.1, 10))
#' 
#' # Call the function to perform cuts
#' newToyDataDF <- performCuts(toyDataDF)
#' 
#' # Clean the environment
#' rm(list=c("toyDataDF", "newToyDataDF"))
#'  
#' @usage performCuts(originalData)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords utilities
#' @export
#
performCuts <- function(originalData) {

  # This function can be tailored to the user needs, if it is necessary to perform
  # data cuts due to missing data bands, or certain flags, etc.
  # originalData <- subset(originalData, originalData$U<30) # Cut points without U magnitude

  return(originalData)
}