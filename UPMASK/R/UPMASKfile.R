#  R package UPMASK file R/UPMASKfile.R
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

#' @title Run UPMASK in a file
#' 
#' @description \code{UPMASKfile} executes the UPMASK method using a file as an input
#' and writes another file as an output. This is a wrapper function that only reads a 
#' file into an R data frame, calls the \code{UPMASKdata} function using this data frame 
#' and the parameters passed by the user and writes the output into another file.
#' 
#' @param filenameWithPathInput a string indicating the file containing the data to run UPMASK on (with full path)
#' @param filenameWithPathOuput a string indicating the file where the output shall be written (with full path) 
#' @param positionDataIndexes an array of integers indicating the columns of the file containing the spatial position measurements
#' @param photometricDataIndexes an array of integers with the column numbers containing photometric measurements (or any other measurement to go into the PCA step)
#' @param photometricErrorDataIndexes an array of integers with the column numbers containing the errors of the photometric measurements
#' @param threshold a double indicating the thresholding level for the random field analysis
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
#' @param fileWithHeader a boolean indicating if the input file has a text header
#' 
#' @references \href{http://dx.doi.org/10.1051/0004-6361/201321143}{Krone-Martins, A. & Moitinho, A., A&A, v.561, p.A57, 2014}
#'
#' @examples
#' \dontrun{
#' # Analyse a simulated open cluster using spatial and photometric data 
#' # Create strings with filenames
#' fileNameI <- "oc_12_500_1000_1.0_p019_0880_1_25km_120nR_withcolors.dat"
#' inputFileName <- system.file("extdata", fileNameI, package="UPMASK")
#' outputFileName <- file.path(tempdir(), "up-RESULTS.dat")
#' 
#' # Example of how to run UPMASK using data from a file
#' # (serious analysis require at least larger nRuns)
#' posIdx <- c(1,2)
#' photIdx <- c(3,5,7,9,11,19,21,23,25,27)
#' photErrIdx <- c(4,6,8,10,12,20,22,24,26,28)
#' UPMASKfile(inputFileName, outputFileName, posIdx, photIdx, photErrIdx, nRuns=5, 
#'            starsPerClust_kmeans=25, verbose=TRUE, fileWithHeader=TRUE)
#' 
#' # Open the resulting file to inspect the results
#' tempResults <- read.table(outputFileName, header=TRUE)
#' 
#' # Create a simple raw plot to see the results
#' pCols <- tempResults[,length(tempResults)]/max(tempResults[,length(tempResults)])
#' plot(tempResults[,1], tempResults[,2], col=rgb(0,0,0,pCols), cex=0.5, pch=19)
#' 
#' # Clean the environment
#' rm(list=c("tempResults", "inputFileName", "outputFileName", "pCols", "fileNameI"))
#' } 
#'  
#' @usage UPMASKfile(filenameWithPathInput, filenameWithPathOuput, 
#' positionDataIndexes=c(1,2), photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
#' photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28), threshold=1, 
#' maxIter=20, starsPerClust_kmeans=50, nstarts_kmeans=50, nRuns=5, 
#' runInParallel=FALSE, paralelization="multicore", independent=TRUE, verbose=FALSE, 
#' autoCalibrated=FALSE, considerErrors=FALSE, finalXYCut=FALSE, fileWithHeader=FALSE)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords misc, utilities
#' @export
#
UPMASKfile <- function(filenameWithPathInput, filenameWithPathOuput, 
					  positionDataIndexes=c(1,2),
					  photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
					  photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28),
					  threshold=1, maxIter=20, starsPerClust_kmeans=50, nstarts_kmeans=50, 
                      nRuns=5, runInParallel=FALSE, paralelization="multicore", 
                      independent=TRUE, verbose=FALSE, autoCalibrated=FALSE, 
                      considerErrors=FALSE, finalXYCut=FALSE, fileWithHeader=FALSE) {
  
  # 
  # User :: Interface
  #
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Starting UPMASK...\n")
    cat(paste("-------------------------------------------------------------------\n"))
  }
  
  # 
  # Data I/O :: Perform File Input
  # 
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Reading the input table from: \n\t",filenameWithPathInput,"\n")
  }
  # Load the file
  ocdata_full <- read.table(filenameWithPathInput, header=fileWithHeader)
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
  }
  
  #
  # Science :: Run the UPMASK caller function
  #
  resultsTable <- UPMASKdata(ocdata_full, 
  					  positionDataIndexes=positionDataIndexes,
					  photometricDataIndexes=photometricDataIndexes,
					  photometricErrorDataIndexes=photometricErrorDataIndexes,
					  threshold=threshold, maxIter=maxIter, 
                      starsPerClust_kmeans=starsPerClust_kmeans, 
                      nstarts_kmeans=nstarts_kmeans, nRuns=nRuns, 
                      runInParallel=runInParallel, paralelization=paralelization, 
                      independent=independent, verbose=verbose, 
                      considerErrors=considerErrors, finalXYCut=finalXYCut)
  
  # 
  # Data I/O :: Perform File Output
  # 
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Writing the output table at: \n\t",filenameWithPathOuput,"\n")
  }
  write.table(resultsTable, filenameWithPathOuput, sep="	", 
              col.names=fileWithHeader, quote=FALSE, row.names=FALSE)	

  # 
  # User :: Interface
  #
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Done!\n")
    cat(paste("-------------------------------------------------------------------\n"))
  }  
}