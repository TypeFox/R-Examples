##
## File:   ProcessXMLFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the XML file processor.
##
processXMLFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("xml"))) return(NULL)

  ## Attempt to read the file.
  xml2::read_xml(path)

}
