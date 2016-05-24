##
## File:   SchemaOnRead.R
## Author: Michael J. North
## Date:   December 22, 2015
##

##
## Define the file existance and extension checker.
##
checkExtensions <- function(path = ".", extensions = NULL) {

  ## Check to make sure that the given path exists.
  if (!file.exists(path)) return(FALSE)

  ## Find the file's extension and normalize it.
  extension <- tolower(tools::file_ext(path))

  ## Check for valid extensions.
  if (is.null(extensions)) return(TRUE)

  ## Check the given file against the requested extensions.
  extension %in% extensions

}
