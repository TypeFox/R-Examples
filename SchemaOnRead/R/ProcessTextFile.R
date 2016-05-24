##
## File:   ProcessTextFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the text file processor.
##
processTextFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("txt"))) return(NULL)

  ## Attempt to read the file.
  read.csv(path, header = FALSE)

}
