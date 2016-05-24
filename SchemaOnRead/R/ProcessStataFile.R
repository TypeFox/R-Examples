##
## File:   ProcessStataFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the Stata file processor.
##
processStataFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("dta"))) return(NULL)

  ## Attempt to read the file.
  haven::read_stata(path)

}
