`ojaSCM` <-
function(X, center= "ojaMedian", p = NULL, silent = FALSE, na.action = na.fail, ...){
    OjaSigns <- ojaSign(X, x = NULL, center = center, p = p, silent = silent, na.action = na.fail, ...)
    SCM <- crossprod(OjaSigns)/nrow(OjaSigns)
    dimnames(SCM) <- NULL
    rownames(SCM) <- colnames(X)
    colnames(SCM) <- colnames(X)
    return(SCM)
}
