lcut3.data.frame <- function(x,
                             context=NULL,
                             name=NULL,
                             parallel=FALSE,
                             ...) {
    if (!is.data.frame(x)) {
        stop("'x' must be a data frame")
    }
    if (ncol(x) <= 0) {
        stop("'x' must contain at least a single column")
    }
    if (!is.null(name)) {
        stop("If 'x' is a matrix or data frame then 'name' must be NULL")
    }
    if (is.null(colnames(x))) {
        stop("Columns of 'x' must have names")
    }
    if (!is.list(context)) {
        context <- rep(list(context), ncol(x))
        names(context) <- colnames(x)
    }
    if (length(intersect(names(context), colnames(x))) != length(names(context))) {
        stop("'context' must be a list with names corresponding to column names")
    }

    loopBody <- function(n) {
        ctx <- context[[n]]
        res <- lcut3(x[, n], 
                    context=ctx,
                    name=n,
                    parallel=FALSE,
                    ...)
        return(res)
    }

    n <- NULL
    if (parallel) {
        result <- foreach(n=colnames(x), .combine=cbind.fsets) %dopar% { return(loopBody(n)) }
    } else {
        result <- foreach(n=colnames(x), .combine=cbind.fsets) %do% { return(loopBody(n)) }
    }
    return(result)
}


lcut5.data.frame <- function(x,
                             context=NULL,
                             name=NULL,
                             parallel=FALSE,
                             ...) {
    if (!is.data.frame(x)) {
        stop("'x' must be a data frame")
    }
    if (ncol(x) <= 0) {
        stop("'x' must contain at least a single column")
    }
    if (!is.null(name)) {
        stop("If 'x' is a matrix or data frame then 'name' must be NULL")
    }
    if (is.null(colnames(x))) {
        stop("Columns of 'x' must have names")
    }
    if (!is.list(context)) {
        context <- rep(list(context), ncol(x))
        names(context) <- colnames(x)
    }
    if (length(intersect(names(context), colnames(x))) != length(names(context))) {
        stop("'context' must be a list with names corresponding to column names")
    }

    loopBody <- function(n) {
        ctx <- context[[n]]
        res <- lcut5(x[, n], 
                    context=ctx,
                    name=n,
                    parallel=FALSE,
                    ...)
        return(res)
    }

    n <- NULL
    if (parallel) {
        result <- foreach(n=colnames(x), .combine=cbind.fsets) %dopar% { return(loopBody(n)) }
    } else {
        result <- foreach(n=colnames(x), .combine=cbind.fsets) %do% { return(loopBody(n)) }
    }
    return(result)
}
