##
## File:   ProcessDBFFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the DBF file processor.
##
processDBFFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("dbf"))) return(NULL)

  ## Attempt to read the file.
  foreign::read.dbf(path)

}
