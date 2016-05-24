#' Utilities for managing a library of SLQ files. 
#' 
#' @name sqlutils-package
#' @aliases sqlutils
#' @docType package
#' @title Utilities for working with SQL files.
#' @author Jason Bryer \email{jason@@bryer.org}
#' @keywords package database sql
#' @import DBI
#' @import roxygen2
#' @import stringr
NULL

#' The locations of SQL files
sqlutils.envir <- new.env()

.onAttach <- function(libname, pkgname) {
	assign("sqlrepos", value=c(paste(system.file(package='sqlutils'), '/sql', sep='')), 
		   envir=sqlutils.envir)
}
