## Return #of iterations for maxim objects

nIter <- function(x, ...)
    ## Number of iterations for iterative models
    UseMethod("nIter")

nIter.default <- function(x, ...)
    x$iterations
