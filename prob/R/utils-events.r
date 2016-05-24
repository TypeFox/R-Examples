

`countrep` <- function (x, ...)
UseMethod("countrep")

`countrep.data.frame` <- function (x, ...){
    apply(x, MARGIN = 1, FUN = countrep, ...)
}


`countrep.default` <- function (x, vals = unique(x), nrep = 2, ...){
    res <- 0
    if (length(x) >= nrep) {
        for (i in 1:length(vals)) {
            if (sum(mapply(all.equal, x, vals[i]) == TRUE) == 
                nrep) {
                res <- res + 1
            }
        }
    }
    return(res)
}


`isin` <- function (x, ...)
UseMethod("isin")



`isin.data.frame` <- function (x, ...){
    apply(x, MARGIN = 1, FUN = isin, ...)
}


`isin.default` <- function (x, y, ordered = FALSE, ...){
    res <- (length(y) <= length(x))
    if (res) {
        temp <- x
        for (i in 1:length(y)) {
            if (is.element(y[i], temp)) {
                if (!ordered) {
                  temp <- temp[-which(temp %in% y[i])[1]]
                }
                else {
                  temp <- temp[-(1:which(temp %in% y[i])[1])]
                }
            }
            else {
                res <- FALSE
                i <- length(y)
            }
        }
    }
    return(res)
}



`isrep` <- function (x, ...)
UseMethod("isrep")


`isrep.data.frame` <- function (x, ...){
    apply(x, MARGIN = 1, FUN = isrep, ...)
}



`isrep.default` <- function (x, vals = unique(x), nrep = 2, ...){
    res <- FALSE
    if (length(x) >= nrep) {
        for (i in 1:length(vals)) {
            if (sum(mapply(all.equal, x, vals[i]) == TRUE) == 
                nrep) {
                res <- TRUE
                i <- length(vals)
            }
        }
    }
    return(res)
}
