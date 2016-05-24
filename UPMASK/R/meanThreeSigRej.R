#  R package UPMASK file R/meanThreeSigRej.R
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
#' @description \code{meanThreeSigRej} will perform an interative rejection using the mean
#' and three sigma of the data. The function will compute the mean and standard 
#' deviation of the input vector, reject all entries lying farther than three-sigma, and 
#' iterate until the mean varies (fractionaly) by less than the tolerance value. The 
#' function will return a data frame with the mean, standard deviation value and the 
#' number of iterations until the convergence was reached.
#' 
#' @param vec a vector with the data to consider
#' @param maxI an integer with the maximum amount of iterations allowed
#' @param tolerance a double with the tolerance value (as \code{(old-new)/old})
#' 
#' @return A data frame with the fields \code{mean}, \code{sd} (the standard deviation) and \code{convergenceAtIter} (the iteration where the convergence was reached).
#' 
#' @examples
#' # Create a simple data set with the values and errors
#' toyData <- c(rnorm(30, mean=0, sd=5), 10000, 1000)
#' 
#' # Call the function to perform cuts
#' toyDataItStat <- meanThreeSigRej(toyData)
#' 
#' cat(paste("True mean             = 0\n"))
#' cat(paste("Before rejection mean =",round(mean(toyData),2),"\n"))
#' cat(paste("After rejetion mean   =",round(toyDataItStat$mean,2),"\n"))
#'
#' # Clean the environment
#' rm(list=c("toyData", "toyDataItStat"))
#'  
#' @usage meanThreeSigRej(vec, maxI, tolerance)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords utilities
#' @export
#
meanThreeSigRej <- function(vec, maxI=50, tolerance=0.001) {
  tvec <- vec
  newM <- mean(tvec)
  i<-0
  for(i in 1:maxI) {
    oldM <- newM
    tvec <- subset(tvec, tvec<(mean(tvec)+3*sd(tvec)))
    newM <- mean(tvec)
    newSd <- sd(tvec)
    # if the mean converged, stop the torture, please...
    if( abs(oldM-newM)/oldM < tolerance) { break }
  }
  if(i==maxI) {
    stop(" ARGH! The Iterative 3-sigma rejection reached the maximum iterations without converging...\n Aborting cowardly!")
  }
  return(data.frame(mean=newM,sd=newSd, convergenceAtIter=i))
}