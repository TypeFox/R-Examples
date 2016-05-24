tail.farules <- function(x, n = 6L, ...) {
    if (!is.farules(x)) {
        stop("'x' is not a valid 'farules' object")
    }
    class(x) <- setdiff(class(x), 'farules')
    return(farules(rules=tail(x$rules, n=n),
                   statistics=tail(x$statistics, n=n)))
}
