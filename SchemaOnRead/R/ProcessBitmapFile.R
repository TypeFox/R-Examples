##
## File:   ProcessBitmapFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the bitmap (BMP, PNG, and JPG) file processor.
##
processBitmapFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("bmp", "png", "jpg"))) return(NULL)

  ## Attempt to read the file.
  readbitmap::read.bitmap(path)

}
