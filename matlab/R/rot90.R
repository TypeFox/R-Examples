###
### $Id: rot90.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Rotates matrix counterclockwise k*90 degrees.
###


##-----------------------------------------------------------------------------
rot90 <- function(A, k = 1) {
    if (!is.matrix(A)) {
        stop(sprintf("argument %s must be matrix", sQuote("A")))
    }

    if (!is.numeric(k)) {
        stop(sprintf("argument %s must be numeric", sQuote("k")))
    } else if (!(length(k) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("k")))
    }

    rot90 <- function(A) {
        n <- matlab::size(A)[2]

        A <- t(A)
        return(A[n:1,])
    }

    rot180 <- function(A) {
        sz <- matlab::size(A)
        m <- sz[1]
        n <- sz[2]

        return(A[m:1, n:1])
    }

    rot270 <- function(A) {
        m <- matlab::size(A)[1]

        return(t(A[m:1,]))
    }

    k <- matlab::rem(k, 4)
    if (k <= 0) {
        k <- k + 4
    }

    return(switch(EXPR = k,
                  rot90(A),
                  rot180(A),
                  rot270(A),
                  A))
}

