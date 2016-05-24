##
## File:   ProcessPajekFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the Pajek file processor.
##
processPajekFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("net", "paj"))) return(NULL)

  ## Attempt to read the file.
  network::read.paj(path)

}
