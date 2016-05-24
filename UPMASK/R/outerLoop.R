#  R package UPMASK file R/outerLoop.R
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

#' @title UPMASK outer loop
#' 
#' @description \code{outerLoop} executes the UPMASK method's outer loop on a data frame, 
#' and returns another data frame as output, with the id of the object and it's 
#' classification as a stellar cluster member or not.
#' 
#' The \code{outerLoop} perform cuts in the data if necessary (by calling 
#' \code{\link{performCuts}}), take errors in the data table into account if the user request (by
#' calling \code{\link{innerLoop}}), runs the inner loop (by calling \code{\link{innerLoop}}) until 
#' convergence of the membership list or util the maximum number of iterations is reached.
#' 
#' @param ocdata_full a data frame with the data to perform the analysis
#' @param positionDataIndexes an array of integers indicating the columns of the data frame containing the spatial position measurements
#' @param photometricDataIndexes an array of integers with the column numbers containing photometric measurements (or any other measurement to go into the PCA step)
#' @param photometricErrorDataIndexes an array of integers with the column numbers containing the errors of the photometric measurements
#' @param threshold a double indicating the thresholding level for the random field analysis
#' @param maxIter an integer the maximum amount of iterations of the outer loop before giving up convergence (usually it is not necessary to modify this)
#' @param plotIter a boolean indicating if the user wants to see iteration plots
#' @param verbose a boolean indicating if the output to screen should be verbose
#' @param starsPerClust_kmeans an integer with the average number of stars per k-means cluster
#' @param nstarts_kmeans an integer the amount of random re-initializations of the k-means clustering method (usually it is not necessary to modify this)
#' @param finalXYCut a boolean indicating if a final cut in the XY space should be performed (defaults to FALSE)
#' @param autoCalibrated a boolean indicating if the number of random field realizations for the clustering check in the position space should be autocalibrated (experimental code, defaults to FALSE).
#' @param considerErrors a boolean indicating if the errors should be taken into account
#' @param run an integer greater than zero indicating the run number
#' @param smartTableDB a database connection to the smart look-up table
#' 
#' @return A data frame with the id and class (member / not member) of each object at this run.
#'
#' @references \href{http://dx.doi.org/10.1051/0004-6361/201321143}{Krone-Martins, A. & Moitinho, A., A&A, v.561, p.A57, 2014}
#' 
#' @examples
#' \dontrun{
#' # Perform a one run of the outerLoop using a simulated open cluster with 
#' # spatial and photometric data 
#' # Load the data into a data frame
#' fileName <- "oc_12_500_1000_1.0_p019_0880_1_25km_120nR_withcolors.dat"
#' inputFileName <- system.file("extdata", fileName, package="UPMASK")
#' ocData <- read.table(inputFileName, header=TRUE)
#' 
#' # Create the look up table
#' library(RSQLite)
#' stcon <- create_smartTable()
#' 
#' # Run the outer loop 
#' posIdx <- c(1,2)
#' photIdx <- c(3,5,7,9,11,19,21,23,25,27)
#' photErrIdx <- c(4,6,8,10,12,20,22,24,26,28)
#' outerLoopRes <- outerLoop(ocData, posIdx, photIdx, PhotErrIdx,
#'                           starsPerClust_kmeans=25, verbose=TRUE, smartTableDB=stcon)
#' 
#' # Clean the environment
#' rm(list=c("inputFileName", "ocData", "posIdx", "photIdx", "photErrIdx", 
#'           "outerLoopRes", "fileName"))
#' dbDisconnect(stcon)
#' } 
#'  
#' @usage outerLoop(ocdata_full, positionDataIndexes=c(1,2),
#' photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
#' photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28), threshold=1, maxIter=25, 
#' plotIter=FALSE, verbose=FALSE, starsPerClust_kmeans=50, nstarts_kmeans=50, 
#' finalXYCut=FALSE, autoCalibrated=FALSE, considerErrors=FALSE, run=0, smartTableDB)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords cluster, methods, multivariate, nonparametric
#' @import stats
#' @export
#
outerLoop <- function(ocdata_full, 
					  positionDataIndexes=c(1,2),
					  photometricDataIndexes=c(3,5,7,9,11,19,21,23,25,27),
					  photometricErrorDataIndexes=c(4,6,8,10,12,20,22,24,26,28),
					  threshold=1, maxIter=25, plotIter=FALSE, 
					  verbose=FALSE, starsPerClust_kmeans=50, nstarts_kmeans=50, 
					  finalXYCut=FALSE, 
					  autoCalibrated=FALSE, considerErrors=FALSE, run=0, 
					  smartTableDB) {
  
  autoThreshold <- TRUE
  
  #library(stats)
  
  # Create an internal index
  idx <- which(ocdata_full$field!=0)
  ocdata_full$field[idx] <- rep(1, times=length(idx))
  
  if(verbose) {
    cat(paste("-------------------------------------------------------------------\n"))
    cat(paste(" Points in the original data   :",length(ocdata_full[,1]),"\n"))
  }
  
  # Perform cuts in the data if necessary
  ocdata_full <- performCuts(ocdata_full)
  
  # Take errors in the data table into account if the user request
  if(considerErrors) {
    ocdata_full <- takeErrorsIntoAccount(ocdata_full, photometricDataIndexes, 
                                         photometricErrorDataIndexes)
  }
  
  if(verbose) {
    cat(paste(" Points after the data cuts    :",length(ocdata_full[,1]),"\n"))
    cat(paste("-------------------------------------------------------------------\n"))
    cat(      " Starting iterations...\n")
    cat(paste(" Maximum iterations :",maxIter,"\n"))
    if(plotIter) {
      cat(paste(" You have choosen to see the iteration's plots.\n"))
    } else {
      cat(paste(" You have choosen NOT to see the iteration's plots.\n"))
    }
    cat(paste("-------------------------------------------------------------------\n"))
  }
  
  # Start iterations
  verb <- 1
  if(!verbose) { 
  	verb <- 0 
  }
  ocdata_full <- data.frame(ocdata_full, id=(1:length(ocdata_full$x))) # create an id
  ocdata_out <- ocdata_full # it start as the full list
  ocdata_res <- ocdata_full # it start as the full list
  solConv <- FALSE
  i <- 0
  for(i in 1:maxIter) {
    oldSize <- length(ocdata_res[,1])
    # Select the columns to analyse
    ocdata <- ocdata_res[,photometricDataIndexes]
    
    # Iterate the inner loop
    ocdata_res <- innerLoop(ocdata_res, ocdata, classAlgol="kmeans", 
                            autoThresholdLevel=threshold, iiter=i, plotIter=FALSE, 
                            starsPerClust_kmeans=starsPerClust_kmeans, 
                            nstarts_kmeans=nstarts_kmeans, verbosity=verb, runId=run, 
                            autoCalibrated=autoCalibrated, 
                            positionDataIndexes=positionDataIndexes, smartTableDB=smartTableDB)
    
    # Flag this iteration's cluster stars
    member <- data.frame(m=rep(0,length(ocdata_out[,1])))
    member$m[ocdata_res$id] <- 1
    ocdata_out <- data.frame(ocdata_out, member)
    
    # Check if the solution converged at this iteration
    if(length(ocdata_res[,1])==oldSize) {
      solConv <- TRUE
      break # not an elegant solution
    }
    
    if(length(ocdata_res[,1])==0) {
      solConv <- FALSE
      break # not an elegant solution
    }
  }
  
  # If the solution converged, then let the user know (if verbose)
  # If it did not, flag everything as non member stars
  if (solConv) {
    if(verbose) {
      cat(paste(" Convergence found at iteration,",i,"!\n"))
    }
  } else {
    if(verbose) {
	    cat(" Sorry, the system never converged...\n But I am not aborting!!\n")
	}
    # Flag this iteration's cluster stars as non-members
    member <- data.frame(m=rep(0,length(ocdata_out[,1])))
    ocdata_out <- data.frame(ocdata_out, member)
  }
  
  # If the user wants to perform a final cut in order to get the spatially clustered
  # members only, then do it.
  if(finalXYCut) {
    if(verbose) {
      cat(paste("-------------------------------------------------------------------\n"))
      cat(      " Selecting the stars at the highest X-Y density levels...\n")
    }
    ocdata_out <- getStarsAtHighestDensityRegion(ocdata_out, threshold=3, verbose=FALSE)
  } else {
    ocdata_out <- data.frame(ocdata_out, finalClass=ocdata_out[,length(ocdata_out)] )	
  }
    
  # That's all folks!  
  return(data.frame(id=ocdata_out$id, class=ocdata_out$finalClass))
}