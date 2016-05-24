##
## File:   ProcessODSFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the OpenDocument spreadsheet file processor.
##
processODSFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("ods"))) return(NULL)

  ## Attempt to read the file.
  readODS::read.ods(file = path)

}
