###
### $Id: isempty.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Determine if object is empty.
###


##-----------------------------------------------------------------------------
isempty <- function(A) {
     return(any(matlab::size(A) == 0))
}

