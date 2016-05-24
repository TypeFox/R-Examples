#  R package UPMASK file R/analyse_randomKde2d_smart.R
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

#' @title Perform analysis of random 2d distributions
#' 
#' @description \code{analyse_randomKde2d_smart} will compute statistics from uniformly 
#' randomly created 2D fields based on Kernel Density Estimations (calling the code 
#' \code{\link{analyse_randomKde2d}}). However, if a random field using the same number of stars
#' was already computed in this run of UPMASK, it will avoid computing it again and will
#' return the value that is stored in a SQLite database table. If the random field was 
#' not yet analysed, it will run the analysis, store the result in the database table, and 
#' return the value.
#' 
#' @param nfields an integer with the number of individual field realisations
#' @param nstars an integer with the number of stars to consider
#' @param maxX the length of the field in X
#' @param maxY the length of the field in Y
#' @param nKde the number of samplings of the kernel in each direction
#' @param showStats a boolean indicating if the user wants to see statistics
#' @param returnStats a boolean indicating if the user wants statistics to be returned
#' @param smartTableDB a database connection to the smart look-up table
#' 
#' @return A data frame with the \code{mean} and \code{sd} fields containing the results 
#' of the random field analysis.
#' 
#' @examples
#' # Create the smart look-up table
#' library(RSQLite)
#' stcon <- create_smartTable()
#' 
#' # Runs the analysis on random fields
#' system.time(
#' toyRes1 <- analyse_randomKde2d_smart(300, 200, 100, 100, smartTableDB=stcon)) # slow
#' system.time(
#' toyRes2 <- analyse_randomKde2d_smart(300, 200, 100, 100, smartTableDB=stcon)) # quick
#' 
#' # Clean the environment
#' rm(list=c("toyRes1", "toyRes2"))
#' dbDisconnect(stcon)
#'  
#' @usage analyse_randomKde2d_smart(nfields=100, nstars, maxX, maxY, nKde=50, 
#' showStats=FALSE, returnStats=TRUE, smartTableDB)
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords utilities
#' @import RSQLite
#' @export
#
analyse_randomKde2d_smart <- function(nfields=100, nstars, maxX, maxY, nKde=50, 
                                      showStats=FALSE, returnStats=TRUE, smartTableDB) {
  
  # Get the table from the database
  smartTable <- dbReadTable(smartTableDB, "smartTable")
  
  res <- subset(smartTable, smartTable$nstar==nstars) # only nstars matters by now...
  
  if(length(res$mean)!=0) {
    retStat <- data.frame(mean=mean(res$mean), sd=mean(res$sd))
  } else {
    # run the analysis...
    retStat <- analyse_randomKde2d(nfields, nstars, maxX, maxY, nKde, showStats, returnStats)
    # and store the results in the table, so next time you won't need to compute it again!
    dbWriteTable(smartTableDB, "smartTable", data.frame(nstar=nstars, mean=retStat$mean, sd=retStat$sd), append=TRUE, row.names = F)
  }
  
  return(retStat)
}