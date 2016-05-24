#  R package UPMASK file R/kde2dForSubset.R
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

#' @title Compute the density based distance quantity using a 2D Kernel Density Estimation
#' 
#' @description \code{kde2dForSubset} will compute the 2D Kernel Density Estimation for 
#' the requested subset of data and will return the quantiy \code{(max(d)-mean(d))/sd(d)}
#' if the option \code{returnDistance} is set to TRUE.
#' 
#' @param df a data frame to use
#' @param setw an integer with the class of the stars to perform the analysis
#' @param n the number of points in the regular grid of the density estimation
#' @param showStats a boolean indicating if the user wants to see output statistics
#' @param printPlots a boolean indicating if the user wants to see plots
#' @param returnDistance a boolean indicating if the distance between the max and the mean in units of standard deviations should be returned
#' @param positionDataIndexes an array of integers indicating the columns of the file containing the spatial position measurements
#'
#' @return A double representing the density based distance quantity.
#'
#' @references \href{http://dx.doi.org/10.1051/0004-6361/201321143}{Krone-Martins, A. & Moitinho, A., A&A, v.561, p.A57, 2014}
#' 
#' @examples
#' # Create a simple data set with the values and errors
#' toyDataDF <- data.frame(x=runif(50, 0, 10), y=runif(50, 0, 10), resMclust.class=rep(1, 50))
#' 
#' # Run the KDE 2D analysis for the required subset
#' disV <- kde2dForSubset(toyDataDF, showStats=FALSE, printPlots=FALSE, returnDistance=TRUE)
#' 
#' # Clean the environment
#' rm(list=c("toyDataDF", "disV"))
#'  
#' @usage kde2dForSubset(df, setw=1, n=50, showStats=TRUE, printPlots=TRUE, 
#' returnDistance=FALSE, positionDataIndexes=c(1,2))
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords utilities
#' @import MASS
#' @export
#
kde2dForSubset <- function(df, setw=1, n=50, showStats=TRUE, printPlots=TRUE, returnDistance=FALSE, positionDataIndexes=c(1,2)) {
  # Load the 2d-kde library
  # library(MASS)
  
  # Filter the data
  dfn <- subset(df, df$resMclust.class==setw)
  
  # Create the 2d-kde
  dataX <- dfn[,positionDataIndexes[1]]
  dataY <- dfn[,positionDataIndexes[2]]
  # using the method of Sheather & Jones (1991) to select the bandwidth
  #	kde2dmap <- kde2d(dataX, dataY, n=n, lims=c(range(df$x), range(df$y)), h = c(width.SJ(dataX), width.SJ(dataY)))
  # using normal bandwidth selection rule of thumb (MASS)
  kde2dmap <- kde2d(dataX, dataY, n=n, lims=c(range(df[,positionDataIndexes[1]]), range(df[,positionDataIndexes[2]])) ) 
  
  # Plot the results
  if(printPlots) {
    image(kde2dmap, main=paste("Class",setw) )
    points(dataX, dataY, pch=19, cex=0.3)
  }
  
  # Print some statistics
  if(showStats) {
    cat(paste("------ Stats - Class",setw," -----\n"))
    cat(paste("   Max dens.   : ", max(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Min dens.   : ", min(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Mean dens.  : ", mean(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Sd dens.    : ", sd(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Median dens.: ", median(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   MAD dens.   : ", mad(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Dist max from the mean (in sd): ", ((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))),"\n" ))
    cat(paste("-------------------------------------\n\n"))		
  }
  
  # Return the distance
  if(returnDistance) {
    return(((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))))
  }
}