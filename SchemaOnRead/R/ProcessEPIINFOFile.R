##
## File:   ProcessEPIINFOFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the EPIINFO file processor.
##
processEPIINFOFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("rec"))) return(NULL)

  ## Attempt to read the file.
  suppressWarnings(foreign::read.epiinfo(path))

}
