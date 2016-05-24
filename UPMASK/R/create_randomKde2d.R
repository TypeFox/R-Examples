#  R package UPMASK file R/create_randomKde2d.R
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
#' @description \code{create_randomKde2d} will compute the 2D Kernel Density Estimation for 
#' a random sampling of the requested number of points and will return the quantiy 
#' \code{(max(d)-mean(d))/sd(d)}, if the option \code{returnDistance} is set to TRUE.
#' 
#' @param nstars an integer with the number of stars to consider
#' @param maxX the length of the field in X
#' @param maxY the length of the field in Y
#' @param nKde the number of samplings of the kernel in each direction
#' @param printPlots a boolean indicating if the user wants to see plots
#' @param showStats a boolean indicating if the user wants to see statistics
#' @param returnDistance a boolean indicating if the user wants statistics to be returned
#' 
#' @return A double representing the density based distance quantity.
#'
#' @references \href{http://dx.doi.org/10.1051/0004-6361/201321143}{Krone-Martins, A. & Moitinho, A., A&A, v.561, p.A57, 2014}
#' 
#' @examples
#' # Compute the density based distance quantity with random fields
#' distVar <- create_randomKde2d(100, 10, 10, showStats=FALSE, 
#'                               printPlots=FALSE, returnDistance=TRUE)
#' 
#' # Clean the environment
#' rm(distVar)
#'  
#' @usage create_randomKde2d(nstars, maxX, maxY, nKde=50, printPlots=FALSE, 
#' showStats=FALSE, returnDistance=FALSE)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords utilities
#' @import MASS
#' @export
#
create_randomKde2d <- function(nstars, maxX, maxY, nKde=50, printPlots=FALSE, showStats=FALSE, returnDistance=FALSE) {
  # Load the 2d-kde library
  #library(MASS)
  
  # Create a random sample in the X-Y space
  dataX <- runif(nstars, 0, maxX)
  dataY <- runif(nstars, 0, maxY)
  
  # Create the 2d-kde
  # using the method of Sheather & Jones (1991) to select the bandwidth
  #	kde2dmap <- kde2d(dataX, dataY, n=nKde, lims=c(range(df$x), range(df$y)), h = c(width.SJ(dataX), width.SJ(dataY)))
  # using normal bandwidth selection rule of thumb (MASS)
  kde2dmap <- kde2d(dataX, dataY, n=nKde, lims=c(0, maxX, 0, maxY) ) 
  
  # Plot the results
  if(printPlots) {
    image(kde2dmap, main=paste("Random Field") )
    points(dataX, dataY, pch=19, cex=0.3)
  }
  
  # Print some statistics
  if(showStats) {
    cat(paste("------ Statistics of the random field -----\n"))
    cat(paste("   Max dens.   : ", max(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Min dens.   : ", min(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Mean dens.  : ", mean(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Sd dens.    : ", sd(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Median dens.: ", median(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   MAD dens.   : ", mad(as.vector(kde2dmap$z)),"\n"))
    cat(paste("   Dist max from the mean (in sd): ", ((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))),"\n" ))
    cat(paste("-------------------------------------------\n\n"))		
  }
  
  # Return the distance
  if(returnDistance) {
    return(((max(as.vector(kde2dmap$z))-mean(as.vector(kde2dmap$z)))/sd(as.vector(kde2dmap$z))))
  }
}