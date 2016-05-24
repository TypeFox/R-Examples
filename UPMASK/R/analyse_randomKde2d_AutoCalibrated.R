#  R package UPMASK file R/analyse_randomKde2d_AutoCalibrated.R
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

#' @title Perform analysis of random 2d distributions (auto calibrated)
#' 
#' @description The experimental code \code{analyse_randomKde2d_AutoCalibrated} 
#' computes statistics from uniformly randomly created 2D fields based on kernel density 
#' estimations (using \code{\link{create_randomKde2d}}). It runs for as many 
#' interations as necessary for the statistical measuremnt stabilize (but it will abort
#' after \code{maxIter} iterations is reached).
#' 
#' @param nstars an integer with the number of stars to consider
#' @param maxX the length of the field in X
#' @param maxY the length of the field in Y
#' @param nKde the number of samplings of the kernel in each direction
#' @param maxIter an integer with the maximum number of iterations
#' @param showStats a boolean indicating if the user wants to see statistics
#' @param returnStats a boolean indicating if the user wants statistics to be returned
#' 
#' @return A data frame with the \code{mean} and \code{sd} fields containing the results 
#' of the random field analysis.
#' 
#' @examples
#' # Runs the analysis on random fields
#' toyRes <- analyse_randomKde2d_AutoCalibrated(100, 100, 100, 100, showStats=TRUE)
#' 
#' # Clean the environment
#' rm(toyRes)
#'  
#' @usage analyse_randomKde2d_AutoCalibrated(nstars, maxX, maxY, nKde=50, maxIter=40, 
#'                                           showStats=FALSE, returnStats=TRUE)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords utilities
#' @import stats
#' @export
#
analyse_randomKde2d_AutoCalibrated <- function(nstars, maxX, maxY, nKde=50, maxIter=40, showStats=FALSE, returnStats=TRUE) {
  
  maxDistStats <- c()
  nTr <- 100
  maxDistTemp <- vector("double", nTr)
  
  iiter <- 0
  
  # run the analysis
  dMaxDist <- 1
  MaxDistOLD <- 0
  while((dMaxDist)>=0.001) {
    for(i in 1:nTr) {
      maxDistTemp[i] <- create_randomKde2d(nstars, maxX, maxY, nKde=nKde, returnDistance=TRUE)
    }
    maxDistStats <- c(maxDistStats, maxDistTemp)
    MaxDisNEW <- mean(maxDistStats)
    dMaxDist <- abs(MaxDisNEW - MaxDistOLD)/MaxDisNEW
    MaxDistOLD <- MaxDisNEW
    iiter <- iiter + 1
    if(iiter>=maxIter) { 
      break
    } 
  }
  
  #	maxDistStats <- vector("double", nfields)
  #	for(i in 1:nfields) {
  #		maxDistStats[i] <- create_randomKde2d(nstars, maxX, maxY, nKde=nKde, returnDistance=TRUE)
  #	}
  
  # Print some statistics
  if(showStats) {
    cat(paste("------ Statistics of the sample random fields -----\n"))
    cat(paste("   Max distance    : ", max(maxDistStats),"\n"))
    cat(paste("   Min distance    : ", min(maxDistStats),"\n"))
    cat(paste("   Mean distance   : ", mean(maxDistStats),"\n"))
    cat(paste("   Sd distance     : ", sd(maxDistStats),"\n"))
    cat(paste("   Median distance : ", median(maxDistStats),"\n"))
    cat(paste("   MAD distance    : ", mad(maxDistStats),"\n"))
    cat(paste("   Iterations      : ", iiter,"\n"))
    cat(paste("---------------------------------------------------\n\n"))		
    hist(maxDistStats, freq=FALSE)
    lines(density(maxDistStats), col="red")
  }
  
  # Return the mean and st.dev of the distribution	
  if(returnStats) {
    return(data.frame(mean=mean(maxDistStats), sd=sd(maxDistStats)))
  }
  
}