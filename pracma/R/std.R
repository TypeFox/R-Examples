##
##  s t d . R
##


std <- function(x, flag=0) {
    if (length(x) == 0) return(c())
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric vector or matrix.")

    n <- if (flag == 0) length(x) - 1 else length(x)
    sqrt(sum((x-mean(x))*(x-mean(x)))/n)
}


std_err <- function(x) {
    if (length(x) == 0) return(c())
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric vector or matrix.")

    sqrt(var(x)/length(x))
}
