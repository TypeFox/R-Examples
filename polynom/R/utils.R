## accumulate a la Abelson and Sussman.
accumulate <- function(f, init, x, right = TRUE) {
    if(length(x) == 0)
        return(init)
    f <- match.fun(f)
    if(right)
        f(x[[1]], Recall(f, init, x[-1], right = TRUE))
    else
        Recall(f, f(init, x[[1]]), x[-1], right = FALSE)
}
