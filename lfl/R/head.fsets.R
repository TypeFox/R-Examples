head.fsets <- function(x, n = 6L, ...) {
    if (!is.fsets(x)) {
        stop("'x' is not a valid 'fsets' object")
    }
    v <- vars(x)
    s <- specs(x)
    class(x) <- setdiff(class(x), 'fsets')
    return(fsets(head(x, n=n), vars=v, specs=s))
}
