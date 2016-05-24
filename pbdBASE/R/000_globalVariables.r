### Do not delete this file and file name!!
### This file should be loaded before all other *.r files.

### This is to avoid the fale positive messages from R CMD check.
###   "no visible binding for global variable"
### Suggested by Prof Brian Ripley
### ?globalVariables

utils::globalVariables(c(".conflicts.OK", ".pbdBASEEnv"))

#' Global Environment for the pbdBASE Package
#' 
#' The environment for the pbdBASE package where "global" variables are stored.
#' 
#' The \code{.__blacs_gridinfo_} and \code{._blacs_initialized} objects are
#' stored in this environment.
#' 
#' @name BASE Global Environment
#' @export
.pbdBASEEnv <- new.env()

.conflicts.OK <- TRUE
