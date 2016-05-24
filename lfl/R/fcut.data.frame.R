fcut.data.frame <- function(x,
                            breaks,
                            name=NULL,
                            type=c('triangle', 'raisedcos'),
                            merge=1,
                            parallel=FALSE,
                            ...) {
    if (!is.data.frame(x)) {
        stop("'x' must be a data.frame")
    }
    if (!is.null(name)) {
        stop("If 'x' is a matrix or data frame then 'name' must be NULL")
    }
    if (is.null(colnames(x))) {
        stop("Columns of 'x' must have names")
    }
    if (!is.list(breaks)) {
        breaks <- rep(list(breaks), ncol(x))
        names(breaks) <- colnames(x)
    }
    if (!is.list(merge)) {
        merge <- rep(list(merge), ncol(x))
        names(merge) <- colnames(x)
    }
    if (!is.list(type)) {
        type <- rep(list(type), ncol(x))
        names(type) <- colnames(x)
    }

    loopBody <- function(n) {
        aBreaks <- breaks[[n]]
        aMerge <- merge[[n]]
        aType <- type[[n]]

        res <- fcut(x[, n], 
                    breaks=aBreaks,
                    name=n,
                    type=aType,
                    merge=aMerge,
                    parallel=FALSE,
                    ...)
        return(res)
    }

    n <- NULL
    if (parallel) {
        result <- foreach(n=colnames(x), .combine=cbind.fsets) %dopar% {
            return(loopBody(n))
        }
    } else {
        result <- foreach(n=colnames(x), .combine=cbind.fsets) %do% {
            return(loopBody(n))
        }
    }
    return(result)
}
