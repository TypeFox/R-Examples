##
## File:   ProcessDIFFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the DIF file processor.
##
processDIFFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("dif"))) return(NULL)

    ## Attempt to read the file.
    try(return(utils::read.DIF(path)), silent = TRUE)

    ## Attempt to read the file again.
    try(utils::read.DIF(path, transpose = TRUE))

}
