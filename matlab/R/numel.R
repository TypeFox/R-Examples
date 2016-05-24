###
### $Id: numel.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Provides number of elements.
###


##-----------------------------------------------------------------------------
numel <- function(A, varargin) {
    if (!missing(varargin)) {
        stop("not implemented")         # need example
    }

    return(prod(matlab::size(A)))
}

