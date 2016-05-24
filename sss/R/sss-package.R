#' Tools for importing files in the triple-s (Standard Survey Structure) format.
#'
#' sss is a set of tools to import survey files in the .sss (triple-s) format.  Triple-s is a standard to transfer survey data between applications.
#' 
#' @references http://www.triple-s.org/
#'
#' The most important exported function is:
#' \code{\link{read.sss}}
#'
#' @docType package
#' @name sss-package
#' @import XML
#' @import  methods 
#' @importFrom utils read.csv type.convert
#' @aliases sss sss-package
#' @keywords package
NULL

.onAttach <- function(libname, pkgname){
  msg <- "The sss package is in early stages of development and still considered experimental."
  packageStartupMessage(msg)
}  

