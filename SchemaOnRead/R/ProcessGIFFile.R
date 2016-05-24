##
## File:   ProcessGIFFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the GIF file processor.
##
processGIFFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("gif"))) return(NULL)

  ## Attempt to read the file.
  caTools::read.gif(path)

}
