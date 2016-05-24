##
## File:   SimpleProcessors.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the simple processors list getter.
##
simpleProcessors <- function() {

  ## Define the simple processors list.
  list(
    processDirectory,
    processDefaultFile
  )

}
