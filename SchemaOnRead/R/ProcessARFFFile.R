##
## File:   ProcessARFFFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the ARFF file processor.
##
processARFFFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("arf", "arff"))) return(NULL)

  ## Attempt to read the file.
  foreign::read.arff(path)

}
