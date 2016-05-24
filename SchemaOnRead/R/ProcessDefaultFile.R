##
## File:   ProcessDefaultFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the default file processor.
##
processDefaultFile <- function(path = ".", ...) {

  ## Check the given path.
  if (file.exists(path)) return(paste(path, "File Type Unknown"))

  ## Return the default value.
  NULL

}
