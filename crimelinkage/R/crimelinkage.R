#' \code{crimelinkage} package: Statistical Methods for Crime Series Linkage
#'
#' Code for criminal case linkage, crime series identification, crime series
#' clustering, and suspect identification.
#'
#' @details The basic inputs will be a data.frame of crime incidents and an offenderTable
#'   data.frame that links offenders to (solved) crimes. 
#' 
#'   The crime incident data must have one column
#'   named \code{crimeID} that provides a unique crime identifier. Other 
#'   recognized columns include: spatial information: \code{X}, \code{Y} which can
#'   be in metric or long/lat; \code{DT.FROM}, \code{DT.TO} for the event times
#'   (these must be of class POSIXct). Other columns containing information about
#'   the crime, crime scene, or suspect can be included as well.
#'   
#'   The offenderTable must have columns: \code{crimeID} (unique crime identifier)
#'   and \code{offenderID} (unique offender identifier).
#'   
#'   See the vignettes for more details.
#'
#' @docType package
#' @name crimelinkage-package
#' @aliases crimelinkage
#' @import graphics stats utils
NULL


