popelem.wprof <-
function(prof, y, ...) {
    if (is.null(dim(y))) y <- t(y)
    res <- apply(y, 1, function(x) {
        tmp <- apply(prof$profiles, 1, function(p) all(p == x))
        if (sum(tmp)==1)
            return(which(tmp))
        else
            return(NA)
    })
    names(res) <- rownames(y)
    return(res)
}
