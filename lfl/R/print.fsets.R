print.fsets <- function(x, ...) {
    if (!is.fsets(x)) {
        stop("'x' is not a valid 'fsets' object")
    }
    v <- vars(x)
    s <- specs(x)
    class(x) <- setdiff(class(x), 'fsets')
    attr(x, 'vars') <- NULL
    attr(x, 'specs') <- NULL
    print(x)
    cat("\nvars:\n")
    print(v)
    cat("\nspecs:\n")
    print(s)
}
