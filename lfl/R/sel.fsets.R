sel.fsets <- function(x, i=rep(TRUE, nrow(x)), j=rep(TRUE, ncol(x)), ...) {
    if (!is.fsets(x)) {
        stop("'x' must be an instance of the 'fsets' class")
    }
    v <- vars(x)
    s <- specs(x)
    class(x) <- setdiff(class(x), 'fsets')
    return(fsets(x[i, j, drop=FALSE], v[j], s[j, j, drop=FALSE]))
}
