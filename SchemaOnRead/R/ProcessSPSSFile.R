##
## File:   ProcessSPSSFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the SPSS file processor.
##
processSPSSFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("sav"))) return(NULL)

  ## Attempt to read the file.
  haven::read_spss(path)

}
