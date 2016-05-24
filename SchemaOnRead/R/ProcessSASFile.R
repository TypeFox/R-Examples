##
## File:   ProcessSASFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Define the SAS file processor.
##
processSASFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("sas"))) return(NULL)

  ## Attempt to read the file.
  haven::read_sas(path)

}
