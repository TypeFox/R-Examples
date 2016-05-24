###
### $Id: mod.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Modulus after division.
###


##-----------------------------------------------------------------------------
mod <- function(x, y) {
    ans <- x %% y
    ## Substitute x[off] in answer anywhere y[off] equals zero
    if (length(zero.off <- which(y == 0))) {
        ans[zero.off] <- if (length(x) == 1) x else x[zero.off]
    }
    ans
}

