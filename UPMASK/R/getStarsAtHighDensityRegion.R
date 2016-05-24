#  R package UPMASK file R/getStarsAtHighestDensityRegion.R
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

#' @title Perform cut in the membership list based on the 2D space distribution
#' 
#' @description \code{getStarsAtHighestDensityRegion} will compute the 2D Kernel Density 
#' Estimation for the requested subset of data and will return only the stars in the most
#' dense region.
#' 
#' @param ocdata_out a data frame to use
#' @param threshold a double with the thresholding level
#' @param posIdx an array of integers indicating the columns of the data frame containing the spatial positions
#' @param plotAnalysis a boolean indicating if the anaylsis should be plotted
#' @param verbose a boolean indicating if the code should be verbose
#' 
#' @return A data frame with the objects which were selected from \code{ocdata_out}
#' 
#' @examples
#' # Create a simple data set
#' toyDataDF <- data.frame(x=runif(50, 0, 10), y=runif(50, 0, 10), resMclust.class=rep(1, 50))
#' toyDataDF <- rbind(toyDataDF, data.frame(x=rnorm(50, 2, 3), 
#'                    y=rnorm(50, 4, 3), resMclust.class=rep(1, 50)))
#' 
#' # Perform the XY density based cut
#' toyRes <- getStarsAtHighestDensityRegion(toyDataDF)
#' 
#' # Clean the environment
#' rm(list=c("toyDataDF", "toyRes"))
#'  
#' @usage getStarsAtHighestDensityRegion(ocdata_out, threshold=2, posIdx=c(1,2), 
#' plotAnalysis=FALSE, verbose=FALSE)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords utilities
#' @import MASS
#' @export
#
getStarsAtHighestDensityRegion <- function(ocdata_out, threshold=2, posIdx=c(1,2), plotAnalysis=FALSE, verbose=FALSE) {
  # Get the cluster stars
  oc_data_clean <- subset(ocdata_out, ocdata_out$resMclust.class==1)
  
  # Perform the Kde2d in the X-Y space
  dens_map <- kde2d(oc_data_clean[,posIdx[1]],oc_data_clean[,posIdx[2]],n=50, 
                    lims=c(range(ocdata_out[,posIdx[1]]), range(ocdata_out[,posIdx[2]])))
  
  # Compute the average density and its sd (using an iterative method)
  dens_vals <- as.vector(dens_map$z)
  stat_dens <- meanThreeSigRej(dens_vals)
  
  # Compute the average density and its sd, difference from the random fields
  #	stat_dens <- analyse_randomKde2d(nfields=2000, nstars=length(oc_data_clean$x), 
  #								(max(oc_data_clean$x)-min(oc_data_clean$x)),
  #								(max(oc_data_clean$y)-min(oc_data_clean$y)), 
  #								nKde=50, showStats=FALSE, returnStats=TRUE)
  
  
  loopFlag <- TRUE
  while(loopFlag) {	
    # Create a density map with flags for the selected cluster / field regions
    if(verbose) {
      print(stat_dens)
    }
    flagged_dens_map <- dens_map
    flagged_dens_map$z[which(flagged_dens_map$z<(stat_dens$mean+threshold*stat_dens$sd))] <- 0
    flagged_dens_map$z[which(flagged_dens_map$z>=(stat_dens$mean+threshold*stat_dens$sd))] <- 1
    
    # Select the stars inside the flagged region
    oc_data_clean_tmp <- oc_data_clean
    # before anything else, reorganize the data
    dmap <- flagged_dens_map
    xvec <- rep(dmap$x,times=length(dmap$y))
    yvec <- vector("double",(length(dmap$x)*length(dmap$y)) )
    kk <- 1
    for (j in 1:length(dmap$y)) {
      yvec[kk:(kk+length(dmap$x)-1)] <- rep(dmap$y[j], length(dmap$x))
      kk <- kk+length(dmap$x)
    }
    zvec <- as.vector(dmap$z)
    flagged_region <- data.frame(x=xvec, y=yvec, z=zvec)
    # first let's get rid of the easy ones...
    flagged_region <- subset(flagged_region, flagged_region$z==1)
    dxP2 <- abs(flagged_dens_map$x[2]-flagged_dens_map$x[1])/2
    dyP2 <- abs(flagged_dens_map$y[2]-flagged_dens_map$y[1])/2
    oc_data_clean_tmp <- subset(oc_data_clean_tmp, oc_data_clean_tmp[,posIdx[1]]>=(min(flagged_region$x)-dxP2))
    oc_data_clean_tmp <- subset(oc_data_clean_tmp, oc_data_clean_tmp[,posIdx[1]]<=(max(flagged_region$x)+dxP2))
    oc_data_clean_tmp <- subset(oc_data_clean_tmp, oc_data_clean_tmp[,posIdx[2]]>=(min(flagged_region$y)-dyP2))
    oc_data_clean_tmp <- subset(oc_data_clean_tmp, oc_data_clean_tmp[,posIdx[2]]<=(max(flagged_region$y)+dyP2))
    
    loopFlag <- FALSE
    
    if(verbose) {
      print(oc_data_clean_tmp)
    }	
    
    if(length(oc_data_clean_tmp[,posIdx[1]])==0) {
      loopFlag=TRUE
      threshold<-threshold-1
      cat(paste(" WARNING: The thresholding for the final cluster selection in the X-Y density map was to high!\n"))
      cat(paste(" WARNING:     since no stars were present in the final cluster, I am lowering the threshold!\n"))
      cat(paste(" WARNING:     The NEW threshold will be",threshold,"\n"))
      if(threshold<0) {
        stop("ABORTING! Threshold for spatial clustering of member stars is less than zero!")
      }
    }
  }
  
  # ok, now check the hard ones one by one
  trueVec <- c()
  for(i in 1:length(oc_data_clean_tmp[,posIdx[1]])) { # for each star
    # check if this star is inside any of the flagged boxes
    isInside <- FALSE
    for(j in 1:length(flagged_region$z)) {
      if( (oc_data_clean_tmp[i,posIdx[1]] >= (flagged_region$x[j]-dxP2) ) &&
            (oc_data_clean_tmp[i,posIdx[1]] <= (flagged_region$x[j]+dxP2) ) &&
            (oc_data_clean_tmp[i,posIdx[2]] >= (flagged_region$y[j]-dyP2) ) &&
            (oc_data_clean_tmp[i,posIdx[2]] <= (flagged_region$y[j]+dyP2) )) {
        isInside <- TRUE
        break 
      }
    }
    if (isInside) {
      trueVec <- c(trueVec, i)
    }
  }
  oc_data_veryclean_tmp <- oc_data_clean_tmp[trueVec,]
  
  # Plot the results
  if(plotAnalysis) {
    clevels <- c((stat_dens$mean), (stat_dens$mean+1*stat_dens$sd), (stat_dens$mean+2*stat_dens$sd), 
                 (stat_dens$mean+3*stat_dens$sd), (stat_dens$mean+threshold*stat_dens$sd))		
    
    # flagged plot
    colvec <- rgb(0:1/2,0,0)
    dev.new()
    plot(ocdata_out[,posIdx[1]], ocdata_out[,posIdx[2]], pch=19, cex=0.1, col="black", type="n", main=paste("Final Open Cluster"), xlab="X", ylab="Y")
    image(flagged_dens_map, col=colvec, add=TRUE)
    contour(dens_map, levels=clevels, labels=c("mean", "1sd", "2sd","3sd", "user"), col=c("red", "red", "red", "red", "blue"), add=TRUE, labcex=1.5)
    points(oc_data_clean[,posIdx[1]], oc_data_clean[,posIdx[2]], pch=19, cex=0.2, col="white") # OC data	
    points(oc_data_clean_tmp[,posIdx[1]], oc_data_clean_tmp[,posIdx[2]], cex=1, col="blue") # OC data	
    points(oc_data_veryclean_tmp[,posIdx[1]], oc_data_veryclean_tmp[,posIdx[2]], pch=19, cex=1, col="green") # OC data	
    
    # normal plot
    collevels <- 256
    div <- 2*collevels
    colvec <- rgb(0:collevels/div,0:collevels/div,0:collevels/div)
    dev.new()
    plot(ocdata_out[,posIdx[1]], ocdata_out[,posIdx[2]], pch=19, cex=0.1, col="black", type="n", main=paste("Final Open Cluster"), xlab="X", ylab="Y")
    image(dens_map, col=colvec, add=TRUE)
    contour(dens_map, levels=clevels, labels=c("mean", "1sd", "2sd","3sd", "user"), col=c("red", "red", "red", "red", "blue"), add=TRUE, labcex=1.5)
    points(oc_data_clean[,posIdx[1]], oc_data_clean[,posIdx[2]], pch=19, cex=0.2, col="white") # OC data
    
  }
  
  retDf <- data.frame(ocdata_out, finalClass=rep(0,length(ocdata_out[,posIdx[1]])))
  
  retDf$finalClass[oc_data_veryclean_tmp$id] <- 1
  
  return(retDf)	
}