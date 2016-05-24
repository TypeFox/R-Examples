#  R package UPMASK file R/UPMASKdata.R
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

#' @title Run UPMASK in a data frame
#' 
#' @description \code{UPMASKdata} executes the UPMASK method on a data frame, and returns 
#' another data frame as output, including the membership analysis result as additional 
#' columns. 
#' 
#' \code{UPMASKdata} is a method for performing membership assignment in stellar 
#' clusters. The distributed code is prepared to use photometry and spatial positions, 
#' but it can take into account other types of data as well. The method is able to take 
#' into account arbitrary error models (the used must rewrite the 
#' \code{\link{takeErrorsIntoAccount}} function), and it is unsupervised, data-driven, 
#' physical-model-free and relies on as few assumptions as possible. The approach followed
#' for membership assessment is based on an iterative process, principal component 
#' analysis, a clustering algorithm and a kernel density estimation.
#' 
#' @param dataTable a data frame with the data to perform the analysis
#' @param positionDataIndexes an array of integers indicating the columns of the data frame containing the spatial position measurements
#' @param photometricDataIndexes an array of integers with the column numbers containing photometric measurements (or any other measurement to go into the PCA step)
#' @param photometricErrorDataIndexes an array of integers with the column numbers containing the errors of the photometric measurements
#' @param threshold a double indicating the thresholding level for the random field analysis
#' @param classAlgol a string indicating the type of clustering algorithm to consider. Only k-means is implemented at this moment (defaults to kmeans)
#' @param maxIter an integer the maximum amount of iterations of the outer loop before giving up convergence (usually it is not necessary to modify this)
#' @param starsPerClust_kmeans an integer with the average number of stars per k-means cluster
#' @param nstarts_kmeans an integer the amount of random re-initializations of the k-means clustering method (usually it is not necessary to modify this)
#' @param nRuns the total number of individual runs to execute the total number of outer loop runs to execute
#' @param runInParallel a boolean indicating if the code should run in parallel
#' @param paralelization a string with the type of paralilization to use. the paralelization can be: "multicore" or "MPIcluster". At this moment only "multicore" is implemented (defaults to multicore).
#' @param independent a boolean indicating if non-parallel runs should be completely independent
#' @param verbose a boolean indicating if the output to screen should be verbose
#' @param autoCalibrated a boolean indicating if the number of random field realizations for the clustering check in the position space should be autocalibrated (experimental code, defaults to FALSE).
#' @param considerErrors a boolean indicating if the errors should be taken into account
#' @param finalXYCut a boolean indicating if a final cut in the XY space should be performed (defaults to FALSE)
#' 
#' @return A data frame with the original data used to run the method and additional columns indicating the classification at each run, as well as a membership probability in the frequentist sense.
#' 
#' @references \href{http://dx.doi.org/10.1051/0004-6361/201321143}{Krone-Martins, A. & Moitinho, A., A&A, v.561, p.A57, 2014}
#'
#' @examples
#' \dontrun{
#' # Analyse a simulated open cluster using spatial and photometric data 
#' # Load the data into a data frame
#' fileNameI <- "oc_12_500_1000_1.0_p019_0880_1_25km_120nR_withcolors.dat"
#' inputFileName <- system.file("extdata", fileNameI, package="UPMASK")
#' ocData <- read.table(inputFileName, header=TRUE)
#' 
#' # Example of how to run UPMASK using data from a data frame
#' # (serious analysis require at least larger nRuns)
#' posIdx <- c(1,2)
#' photIdx <- c(3,5,7,9,11,19,21,23,25,27)
#' photErrIdx <- c(4,6,8,10,12,20,22,24,26,28)
#' 
#' upmaskRes <- UPMASKdata(ocData, posIdx, photIdx, PhotErrIdx, nRuns=2, 
#'                         starsPerClust_kmeans=25, verbose=TRUE)
#' 
#' # Create a simple raw plot to see the results
#' pCols <- upmaskRes[,length(upmaskRes)]/max(upmaskRes[,length(upmaskRes)])
#' plot(upmaskRes[,1], upmaskRes[,2], col=rgb(0,0,0,pCols), cex=0.5, pch=19)
#' 
#' # Clean the environment
#' rm(list=c("inputFileName", "ocData", "posIdx", "photIdx", "photErrIdx", 
#'           "upmaskRes", "pCols"))
#' } 
#'  
#' @usage UPMASKdata(dataTable, positionDataIndexes=c(1,2),
#' photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
#' photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28), threshold=1, 
#' classAlgol="kmeans", maxIter=25, starsPerClust_kmeans=25, nstarts_kmeans=50, 
#' nRuns=8, runInParallel=FALSE, paralelization="multicore", independent=TRUE, 
#' verbose=FALSE, autoCalibrated=FALSE, considerErrors=FALSE, finalXYCut=FALSE)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords cluster, methods, multivariate, nonparametric
#' @import parallel
#' @export
#
UPMASKdata <- function(dataTable, 
					  positionDataIndexes=c(1,2),
					  photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
					  photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28),
					  threshold=1, classAlgol="kmeans", maxIter=25, 
					  starsPerClust_kmeans=25, nstarts_kmeans=50, nRuns=8, 
					  runInParallel=FALSE, paralelization="multicore", independent=TRUE, 
					  verbose=FALSE, autoCalibrated=FALSE, considerErrors=FALSE, 
					  finalXYCut=FALSE) {
  
  # 
  # User :: Interface
  #
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Starting UPMASK analysis!\n")
    cat(      " UPMASK kernels v. 1.0\n")
    cat(paste("-------------------------------------------------------------------\n"))
  }
  
  # 
  # Science
  #
  # Create the smart look-up table for the random analysis
  # Unfortunatelly this must be a global variable so it can be shared among 
  # parallel processes
  stcon <- create_smartTable() 
    
  # Create a place-holder list
  pp <- list()
  
  # Run UPMASK outerloop (and innerloop, called inside the outerloop)
  if (!runInParallel) {
    if(independent) {
      for(i in 1:nRuns) {
        pp[[length(pp)+1]] <- outerLoop(dataTable, threshold=threshold,
        						maxIter=maxIter, plotIter=FALSE, 
        						starsPerClust_kmeans=starsPerClust_kmeans, 
                                nstarts_kmeans=nstarts_kmeans,  
                                verbose=verbose, finalXYCut=finalXYCut, 
                                autoCalibrated=autoCalibrated, 
                                considerErrors=considerErrors, run=i, smartTableDB=stcon)
      }
    }
  } else {
    if(paralelization=="multicore") {
      # If the user wants to run in an SMP machine, then we want to use occupy the cores!
      # library(parallel)
      # Ok, now it is just a matter of running the code using a list
      pp <- mclapply(1:nRuns, function(x) { outerLoop(dataTable, threshold=threshold, 
                                              maxIter=maxIter, plotIter=FALSE, 
                                              starsPerClust_kmeans=starsPerClust_kmeans, 
                                              nstarts_kmeans=nstarts_kmeans, 
                                              verbose=verbose, finalXYCut=finalXYCut, 
                                              autoCalibrated=autoCalibrated, 
                                              considerErrors=considerErrors, run=x, 
                                              smartTableDB=stcon)} )
    } else if(paralelization=="MPIcluster") {
      # if the user wants to run in a cluster
      # then we need to use MPI
      stop(" MPIcluster is not implemented yet.\n Aborting cowardly!")
    }
  }
  
  ## To do :: check if this is necessary, after all the user
  ## should perform the cuts himself
  dataTable <- performCuts(dataTable) 
  
  # Organize the results into a single data frame
  mergedResults <- data.frame(id=1:length(dataTable[,1]))
  for(i in 1:length(pp)) {
    if(length(pp[[i]])>1) {
      mergedResults <- data.frame(mergedResults, class=pp[[i]]$class)
    }
  }
  
  # Compute a frequency to output as a frequentist probability
  freq <- vector("double",length(mergedResults$id))
  nConverged <- length(mergedResults)-1
  if(nConverged > 0) {
    for(i in 1:length(mergedResults$id)) {
      freq[i] <- sum(mergedResults[i,2:(nConverged+1)])/nConverged
    }
    mergedResults <- cbind(mergedResults, probability=freq)
  } else {
    for(i in 1:length(mergedResults$id)) {
      freq[i] <- 0
    }
    mergedResults <- cbind(mergedResults, probability=freq)
  }
  
  # Organize the results into a single data frame
  ocdata_out <- data.frame(dataTable, mergedResults[,2:length(mergedResults)])
  
  # Clean the global variable (argh!) storing the smart lookup table...
  #rm(smartTable, pos=globalenv())
  
  # Close the connection to the database storing the smart lookup table
  dbDisconnect(stcon)
  
  # 
  # User :: Interface
  #
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " UPMASK analysis is finished!\n")
    cat(paste("-------------------------------------------------------------------\n"))
  }
  	
  # That's all folks!
  return(ocdata_out)
}