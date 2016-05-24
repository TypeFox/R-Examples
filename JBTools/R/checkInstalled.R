checkInstalled <- function(
  ##title<< Check whether a command can be invoked via the command line
  commandName ##<< character string: name of the program/command to check
  )
  ##description<<
  ## checkInstalled checks whether an external command can be run on the command line.
  ##details<<
  ## The test is a simple wrapper around Sys.which which returns TRUE if which returns
  ## a character string and FALSE if not. 
{
 out <- Sys.which(commandName)
 if (length(out) == 0) {
   output <- FALSE
 } else {
   output <- TRUE
 }
 ##value<< logical: whether the program is installed.
 return(output)
}
