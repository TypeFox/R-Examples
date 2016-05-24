##
## File:   ProcessSystatFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the Systat file processor.
##
processSystatFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("sys"))) return(NULL)

  ## Attempt to read the file.
  foreign::read.systat(path)

}
