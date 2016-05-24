#' @title Internal functions used in FSAdata.
#' 
#' @description Internal functions used in FSAdata.
#' 
#' @rdname FSAdata-internals
#' 
#' @keywords internal
#' 
#' @aliases .onAttach
#' 
.onAttach <- function(lib,pkg,...) {
  ## Get version number -- basically code from readPkgVersion in SweaveListingUtils
  vers <- read.dcf(system.file("DESCRIPTION",package=pkg,lib.loc=lib),fields="Version")
  vers <- paste0(vers,"         ",ifelse(nchar(vers)==5," ",""))
  ## Send message
  msg <- paste("\n\n")
  msg <- paste(msg,"##############################################\n")
  msg <- paste(msg,"##","FSAdata package, version",vers,"##\n")
  msg <- paste(msg,"##   by Derek H. Ogle, Northland College    ##\n")
  msg <- paste(msg,"##                                          ##\n")
  msg <- paste(msg,"## Run ?FSAdata for documentation with      ##\n")
  msg <- paste(msg,"##   search hints to find data for specific ##\n")
  msg <- paste(msg,"##   types of fisheries analyses.           ##\n")    
  msg <- paste(msg,"##############################################\n\n")
  packageStartupMessage(msg)
}