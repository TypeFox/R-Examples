#  R package UPMASK file R/innerLoop.R
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

#' @title UPMASK inner loop
#' 
#' @description \code{innerLoop} executes the UPMASK method's inner loop and returns 
#' the stars which were considered as cluster member stars.
#' 
#' The \code{innerLoop} perform the PCA, runs the clustering algorithm and check for 
#' overdensities in the spatial distribution for the clustered stars in the PC space using
#' a 2d kernel density estimation. 
#' 
#' @param ocdata_full a data frame with the data to perform the analysis
#' @param ocdata a data frame with the data to consider in the PCA step
#' @param classAlgol a string indicating the type of clustering algorithm to consider. Only k-means is implemented at this moment (defaults to kmeans)
#' @param autoThresholdLevel an integer indicating the level for thresholding of the spatial distribution
#' @param autoThreshold a boolean indicating if autoThresolding should be adopted (defaults to TRUE)
#' @param iiter and integer indicating the number of the iteration (passed by the \code{outerLoop})
#' @param plotIter a boolean indicating if the user wants to see iteration plots (defaults to FALSE)
#' @param verbosity a flag indicating the verbosity level: it can be 0 (no screen output at all), 1 (minimum), >=2 (all)
#' @param starsPerClust_kmeans an integer with the average number of stars per k-means cluster
#' @param nstarts_kmeans an integer the amount of random re-initializations of the k-means clustering method (usually it is not necessary to modify this)
#' @param runId an integer greater than zero indicating the run Id (passed by the \code{outerLoop})
#' @param autoCalibrated a boolean indicating if the number of random field realizations for the clustering check in the position space should be autocalibrated (experimental code, defaults to FALSE).
#' @param stopIfEmpty a boolean indicating if the code should completely stop if no spatial clustering is detected (defaults to FALSE)
#' @param positionDataIndexes an array of integers indicating the columns of the data frame containing the spatial position measurements
#' @param smartTableDB a database connection to the smart look-up table
#'
#' @return A data frame with objects considered as members at this iteration.
#'
#' @references \href{http://dx.doi.org/10.1051/0004-6361/201321143}{Krone-Martins, A. & Moitinho, A., A&A, v.561, p.A57, 2014}
#' 
#' @examples
#' \dontrun{
#' # Perform a one run of the innerLoop using a simulated open cluster with 
#' # spatial and photometric data 
#' # Load the data into a data frame
#' fileName <- "oc_12_500_1000_1.0_p019_0880_1_25km_120nR_withcolors.dat"
#' inputFileName <- system.file("extdata", fileName, package="UPMASK")
#' ocData <- read.table(inputFileName, header=TRUE)
#' ocData <- data.frame(ocData, id=(1:length(ocData[,1]))) # create an id
#' 
#' # Prepare the data to run the inner loop
#' posIdx <- c(1,2)
#' photIdx <- c(3,5,7,9,11,19,21,23,25,27)
#' 
#' # Create the look up table
#' library(RSQLite)
#' stcon <- create_smartTable()
#' 
#' # Run the inner loop 
#' innerLoopRes <- innerLoop(ocData, ocData[,photIdx], autoThresholdLevel=1, verbosity=2,
#'                           starsPerClust_kmeans=25, positionDataIndexes=posIdx, 
#'                           smartTableDB=stcon)
#' 
#' # Clean the environment
#' rm(list=c("inputFileName", "ocData", "posIdx", "photIdx", "innerLoopRes", 
#'    "fileName"))
#' dbDisconnect(stcon)
#' } 
#'  
#' @usage innerLoop(ocdata_full, ocdata, classAlgol="kmeans", autoThresholdLevel=3, 
#' autoThreshold=TRUE, iiter=0, plotIter=FALSE, verbosity=1, starsPerClust_kmeans=50, 
#' nstarts_kmeans=50, runId=0, autoCalibrated=FALSE, stopIfEmpty=FALSE, 
#' positionDataIndexes=c(1,2), smartTableDB)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords cluster, methods, multivariate, nonparametric
#' @export
#
innerLoop <- function(ocdata_full, ocdata, classAlgol="kmeans", autoThresholdLevel=3, 
                         autoThreshold=TRUE, iiter=0, plotIter=FALSE, verbosity=1, 
                         starsPerClust_kmeans=50, nstarts_kmeans=50, runId=0, 
                         autoCalibrated=FALSE, stopIfEmpty=FALSE, 
                         positionDataIndexes=c(1,2), smartTableDB) {
  if(verbosity!=0) {
    cat(paste(" [runId:",runId,"] ITERATION:",iiter," RUNNING...\n"))
  }
  inSize <- length(ocdata_full)
  
  # Perform a PCA
  if(verbosity>=1) {
    cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- [1/3] Performing PCA...\n"))
  }
  ocdata_pca <- prcomp(ocdata, scale=TRUE, center=TRUE, cor=TRUE)
  ocdata_px <- predict(ocdata_pca)
  
  # Perform the clustering
  if (classAlgol=="kmeans") {
    if(verbosity>=1) {
      cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- [2/3] Performing k-means clustering analysis in the PCA space...\n"))
    }
    starsPerClust <- starsPerClust_kmeans
    nclust <- round(length(ocdata_full[,1])/starsPerClust)
    if(nclust > 1) {
      fit <- kmeans(ocdata_px[,1:4], nclust, nstart=nstarts_kmeans, iter.max=100)
      # get cluster means
      aggregate(ocdata_px, by=list(fit$cluster), FUN=mean)
      # append cluster assignment
      ocdata_px <- data.frame(ocdata_px, resMclust.class=fit$cluster)
    } else {
      # if the number of predicted clusters is less than one, lets prevent the kmeans method from 
      # crashing the code by assigning a randomized choice as the result (this is expected to happen if
      # a very small number of stars are assigned at a certain iteration).
      ocdata_px <- data.frame(ocdata_px, resMclust.class=round(runif(length(ocdata_px[,1]), 0, 1)))
    } 
  } else {
     stop(" Error: the selected method for the clustering in the inner loop is not implemented.")
  }
  
  # Print PCA info and create plots -- KEPT FOR DEPURATION PURPOSES
  if(plotIter){
    dev.new()
    par(cex=0.3)
    #	plot(resMclust, data=ocdata_px[,1:4])
    pairs(ocdata_px[,1:4], pch=19, cex=0.2, col=rainbow(max(ocdata_px$resMclust.class))[ocdata_px$resMclust.class])
  }
  
  # Merge the colors, the PCA columns and cluster classification with the original results
  ocdata_full_withColorsAndPcaCols <- data.frame(ocdata_full, ocdata, ocdata_px)
    
  # Select the classes with densities above the threshold
  if(verbosity>=1) {
    cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- [3/3] Performing comparison with random density fields for",max(ocdata_px$resMclust.class),"individual classes...\n"))
    if(verbosity==1) {
      cat(paste(" [runId:",runId,"] ITERATION:",iiter,"          You can get a coffee. This can take long! \n"))
    }
  }
  
  vclass <- c()
  not_class <- c()
  for(i in 1:max(ocdata_px$resMclust.class)) {
    dfn <- subset(ocdata_full_withColorsAndPcaCols, ocdata_full_withColorsAndPcaCols$resMclust.class==i)
        
    if(length(dfn$x)>2) {
      
      dif_max_mean <- kde2dForSubset(ocdata_full_withColorsAndPcaCols, setw=i, returnDistance=TRUE, showStats=FALSE, printPlots=FALSE, positionDataIndexes=positionDataIndexes)
      
      # First, get the thresholding level...
      if(autoThreshold) {
        if(verbosity>=2) {
          cat(paste(" [runId:",runId,"] ITERATION:",iiter,"	-- Class",i," -- Performing analysis of random fields...\n"))
        }
        
        if(autoCalibrated) {
          # This is an experimental code for performing automatic calibration of the 
          # number of random realisations
          at <- analyse_randomKde2d_AutoCalibrated(
                 nstars=length(dfn$resMclust.class), 
                 (max(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[1]])-min(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[1]])),
                 (max(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[2]])-min(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[2]])), 
                 nKde=50, showStats=FALSE, returnStats=TRUE)
        } else {					
          at <- analyse_randomKde2d_smart(
                 nfields=2000, nstars=length(dfn$resMclust.class), 
                 (max(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[1]])-min(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[1]])),
                 (max(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[2]])-min(ocdata_full_withColorsAndPcaCols[,positionDataIndexes[2]])), 
                 nKde=50, showStats=FALSE, returnStats=TRUE, smartTableDB=smartTableDB)
        }
        threshold <- at$mean + autoThresholdLevel*at$sd
        if(verbosity>=2) {
          cat(paste(" [runId:",runId,"] ITERATION:",iiter," 		Automatic threshold for class",i," spatial clustering selected at ",round(threshold,1),"above the mean density.\n"))
          cat(paste(" [runId:",runId,"] ITERATION:",iiter," 		                        class",i," got dif_max_mean             = ",round(dif_max_mean,1),"above the mean density.\n"))
        }
      } else {
      	stop(" The code without autothresolding was deprecated. \b Aborting cowardly.\n\n")
      }
      
      if(is.na(round(threshold,1))) {
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"/- PROBLEM REPORT -----------------------------------------------------------------\n"))
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"| 		Class",i," has a NA value in the threshold!\n"))
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"| 		probably due to its small number of stars: ",length(dfn$x),".\n"))
        kde2dForSubset(ocdata_full_withColorsAndPcaCols, setw=i, returnDistance=FALSE, showStats=TRUE, printPlots=FALSE, positionDataIndexes=positionDataIndexes)
        print(dfn)
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"\\----------------------------------------------------------------------------------\n"))
      } 
      
      if(!is.na(round(threshold,1))) {
        if(round(dif_max_mean,1) >= round(threshold,1)) {
          vclass <- c(vclass, i)
          if(verbosity>=2) {
            cat(paste(" [runId:",runId,"] ITERATION:",iiter," <<<< --	Class",i," : ok!\n"))
          }
        } else {
          not_class <- c(not_class, i)
          if(verbosity>=2) {
            cat(paste(" [runId:",runId,"] ITERATION:",iiter,"      --	Class",i," : WILL BE ELIMINATED!!\n"))
          }
        }
      }
    } else {
      if(verbosity>=2) {
        cat(paste(" [runId:",runId,"] ITERATION:",iiter,"	-- Class",i," -- THERE ARE TWO OR LESS STARS IN THIS CLASS! \n"))
      }
      not_class <- c(not_class, i)
    }
  }
  
  if(length(vclass)==0) {
    cat(" No spatial clustering detected in the real space based on the clustered photometric data in the PCA space!\n")
    #cat(paste(" Spatial thresholding at ",round(threshold,2), "\n")) ### ---- USED FOR DEPURATION PURPOSES
    if(stopIfEmpty) {
      stop(" No spatial clustering detected!\n Aborting cowardly!")
    } else {
      oc_reconst <- ocdata_full_withColorsAndPcaCols[0,]
    }
  }
  
  if(verbosity>=1) {
    if(length(vclass)!=0) {
      if(verbosity >= 2) {
         cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- The ",length(vclass),"selected classes from the kernel density estimation in the X-Y space are: "))
         print(vclass)
      } else {
         cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- Number of classes selected from the kde analysis : ",length(vclass)))
      }
      cat("\n")
    } else {
      cat(paste(" [runId:",runId,"] ITERATION:",iiter," -- No classes were selected from the kernel density estimation in the X-Y space.\n"))
    }
  }
  
  # Organize the data for returning to the outer loop
  if(length(vclass)>=1) {
    # First get the data of the selected objects
    for(i in 1:length(vclass)) {
      oc_tmp <- subset(ocdata_full_withColorsAndPcaCols, ocdata_full_withColorsAndPcaCols$resMclust.class==vclass[i])
      if(i==1) {
        oc_reconst <- oc_tmp
      } else {
        oc_reconst <- rbind(oc_reconst, oc_tmp)
      }
    }
    # now the data of the not selected objects
    for(i in 1:length(not_class)) {
      not_tmp <- subset(ocdata_full_withColorsAndPcaCols, ocdata_full_withColorsAndPcaCols$resMclust.class==not_class[i])
      if(i==1) {
        field_reconst <- not_tmp
      } else {
        field_reconst <- rbind(field_reconst, not_tmp)
      }			
    }
  }
  
  if(verbosity!=0) {
    cat(paste(" [runId:",runId,"] ITERATION:",iiter," DONE!\n"))
    cat(paste("-------------------------------------------------------------------\n"))
  }
  
  # What is going out of this function must be only the astronomical cluster stars only...
  return(oc_reconst[,1:inSize])    
}