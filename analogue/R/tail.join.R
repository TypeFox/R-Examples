## show the last few lines of each joined data set
`tail.join` <- function(x, ...) {
    if(inherits(x, "data.frame"))
        NextMethod("join")
    else
        lapply(x, tail, ...)
}
