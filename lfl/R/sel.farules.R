sel.farules <- function(x, i=rep(TRUE, nrow(x)), ...) {
    if (!is.farules(x)) {
        stop("'x' must be an instance of the 'farules' class")
    }
    return(farules(x$rules[i, drop=FALSE], x$statistics[i, , drop=FALSE]))
}
