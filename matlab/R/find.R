###
### $Id: find.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Find indices of nonzero elements.
###


##-----------------------------------------------------------------------------
find <- function(x) {
    expr <- if (is.logical(x)) {
                x
            } else {
                x != 0
            }
    return(which(expr))
}

