###
### $Id: ndims.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Provides the number of dimensions.
###


##-----------------------------------------------------------------------------
ndims <- function(A) {
    return(length(matlab::size(A)))
}

