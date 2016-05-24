## simple method to head each component of join
`head.join` <- function(x, ...) {
    if(inherits(x, "data.frame"))
        NextMethod("join")
    else
        lapply(x, head, ...)
}
