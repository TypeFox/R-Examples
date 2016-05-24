#  R package UPMASK file R/create_smartTable.R
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

#' @title Create a look up table
#' 
#' @description \code{create_smartTable} will create a look up table for random 
#' field analysis inside an SQLite database. The table is automatically filled each time 
#' \code{UPMASK} calls the function \code{\link{analyse_randomKde2d_smart}}. 
#' 
#' @return A data base connection to the SQLite database containing the smartTable.
#' 
#' @examples
#' # Create the table
#' library(RSQLite)
#' stcon <- create_smartTable() 
#' 
#' # Clean the environment
#' dbDisconnect(stcon)
#'
#' @usage create_smartTable()
#' 
#' @author Alberto Krone-Martins, Andre Moitinho
#' 
#' @keywords utilities
#' @import RSQLite DBI
#' @export
#
create_smartTable <- function() {
#  con <- dbConnect(RSQLite::SQLite(), ":memory:")
  con <- dbConnect(dbDriver("SQLite"), dbname = tempfile())
  dbWriteTable(con, "smartTable", data.frame(nstar=c(0), mean=c(0), sd=c(0)), row.names = F)
  return(con)
}
