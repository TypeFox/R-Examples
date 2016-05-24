##
## File:   ProcessRDSFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the RDS file processor.
##
processRDSFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("rds"))) return(NULL)

  ## Attempt to read the file.
  base::readRDS(path)

}
