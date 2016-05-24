varUsed <- function(x, by.tree=FALSE, count=TRUE) {
    if (!inherits(x, "snpRF"))
        stop(deparse(substitute(x)), "is not a snpRF object")
    if (is.null(x$forest))
        stop(deparse(substitute(x)), "does not contain forest")
	
    ### modified this ###    

    p<- length(x$forest$ncat.autosome) + (length(x$forest$ncat.xchrom)/2) # Total number of variables.

    #####################

    if (count) {
        if (by.tree) {
            v <- apply(x$forest$bestvar, 2, function(x) {
                xx <- numeric(p)
                y <- table(x[x>0])
                xx[as.numeric(names(y))] <- y
                xx
            })
        } else {
            v <- numeric(p)
            vv <- table(x$forest$bestvar[x$forest$bestvar > 0])
            v[as.numeric(names(vv))] <- vv
        }
    } else {
        v <- apply(x$forest$bestvar, 2, function(x) sort(unique(x[x>0])))
        if(!by.tree) v <- sort(unique(unlist(v)))
    }
    v
}
