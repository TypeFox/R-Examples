##
## INPUT:
## b: parameter vector
## regr: list of regressor matrices
## nutils: number of utility equations
## unames: names of utility equations
## finit: whether to coerce utilities not to contain any Infs
##
## RETURN:
## list of numeric vectors (named according to 'unames') of fitted utilities,
## along with element 'b' containing unused parameters (those pertaining to
## variance terms)
## 
makeUtils <- function(b, regr, nutils, unames, finit = TRUE)
{
    utils <- vector("list", nutils)
    if (!missing(unames)) {
        if (length(unames) != nutils)
            stop("length(unames) must equal nutils")
        names(utils) <- unames
    }

    rcols <- sapply(regr, ncol)
    for (i in 1:nutils) {
        if (rcols[i] > 0) {
            utils[[i]] <- as.numeric(regr[[i]] %*% b[1:rcols[i]])
            if (finit)
                utils[[i]] <- finitize(utils[[i]])
            b <- b[-(1:rcols[i])]
        } else {
            utils[[i]] <- rep(0, nrow(regr[[i]]))
        }
    }

    utils$b <- b
    return(utils)
}
