##
## File:   ProcessTIFFFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the TIFF file processor.
##
processTIFFFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("tif", "tiff"))) return(NULL)

  ## Attempt to read the file.
  tiff::readTIFF(path)

}
