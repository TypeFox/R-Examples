###
### $Id: pascal.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Create Pascal matrix.
###


##-----------------------------------------------------------------------------
pascal <- function(n, k = 0) {
    if (!is.numeric(n)) {
        stop(sprintf("argument %s must be numeric", sQuote("n")))
    } else if (!(length(n) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("n")))
    }

    if (!is.numeric(k)) {
        stop(sprintf("argument %s must be numeric", sQuote("k")))
    } else if (!(length(k) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("k")))
    }

    stopifnot(k >= 0, k <= 2)

    P <- diag((-1) ^ seq(0, as.double(n)-1))
    P[,1] <- matlab::ones(n, 1)

    ## Generate Pascal Cholesky factor
    for (j in seq(2, n-1)) {
        for (i in seq(j+1, n)) {
            P[i, j] <- P[i-1, j] - P[i-1, j-1]
        }
    }

    if (k == 0) {
        P <- P %*% t(P)
    } else if (k == 1) {
        ;
    } else if (k == 2) {
        P <- matlab::rot90(P, 3)
        if ((n / 2) == round(n / 2)) {
            P <- -P
        }
    }

    return(P)
}

