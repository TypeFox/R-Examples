##
## File:   ProcessCSVFile.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the CSV file processor.
processCSVFile <- function(path = ".", ...) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("csv"))) return(NULL)

  ## Attempt to read the file.
  utils::read.csv(path, header = TRUE)

}
