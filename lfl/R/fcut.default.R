fcut.default <- function(x, 
                         name=deparse(substitute(x)),
                         ...) {
    if (is.null(x)) {
        return(NULL)
    }
    if (is.null(name)) {
        stop('"name" must not be NULL')
    }

    d <- data.frame(x=x)
    colnames(d) <- paste(name, '.', sep='')
    res <- model.matrix(~ . + 0, data=d)

    # remove unnecessary attributes from the result
    attributes(res)$assign <- NULL
    attributes(res)$contrasts <- NULL
    res <- as.matrix(res)

    theVars <- rep(name, ncol(res))
    names(theVars) <- colnames(res)

    return(fsets(res,
                 vars=theVars,
                 specs=matrix(0,
                              nrow=ncol(res),
                              ncol=ncol(res),
                              dimnames=list(colnames(res), colnames(res)))))
}
