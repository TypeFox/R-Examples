##
## File:   ProcessHTMLFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the HTML file processor.
##
processHTMLFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("htm", "html"))) return(NULL)

  ## Attempt to read the file.
  xml2::read_html(path)

}
