###
### $Id: logspace.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Generate logarithmically spaced vectors.
###


##-----------------------------------------------------------------------------
logspace <- function(a, b, n = 50) {
    if (b == pi) {
        b <- log10(pi)
    }

    10 ^ matlab::linspace(a, b, n)
}

