`mefaNestedless` <-
function(x)
{
    inherits(x, "mefa") || stop("'x' must be of class 'mefa'")
    if (is.null(x$segm))
        stop("no segments in 'x'")
    out <- x$segm
    n <- length(out)
    nam <- nam2 <- names(out)
    for (i in 1:(nchar(nam2[1]) + 1))
        substr(nam2[-1], 1, nchar(nam2[1]) + 1) <- ""
    names(out) <- nam2
    for (i in 2:n) {
        out[[i]] <- out[[i]] - x$segm[[(i-1)]]
    }
    x$segm <- out
    x$call <- match.call()
    attr(x, "nested") <- FALSE
    return(x)
}
