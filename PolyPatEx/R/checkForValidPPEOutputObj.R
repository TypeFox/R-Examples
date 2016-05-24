
## Function checkForValidPPEOutputObj
##
## A utility function to perform a (very simple) check that an object
## passed through is a \code{\link{genotPPE}} or
## \code{\link{phenotPPE}} output object.  If the object fails the
## check, an error is produced via \code{\link{stop}} - otherwise the
## function returns \code{\link{invisible(NULL)}}
##
## @title Check for a genotPPE/phenotPPE output object
## @param x an object to be checked.
## @return If no error is produced, this function returns NULL (invisibly)
## @author Alexander Zwart (alec.zwart at csiro.au)
## @examples
## \dontrun{
##
## checkForValidPPEDataset(x)
##
## }
##
checkForValidPPEOutputObj <- function(x)
  {
    if (!all(names(x) %in% c("progenyTables","adultTables")))
      {
        errorString <- paste("\n\n This function should be passed the object returned by function",
                             "\n genotPPE or function phenotPPE().  You do not appear to have ",
                             "\n passed through the correct type of object...\n\n",
                             collapse=" ")
        stop(errorString)
      } else {
        return(invisible(NULL))
      }
  }
