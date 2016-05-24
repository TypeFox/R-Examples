typeII2typeI <- function(x)
{
    if (!inherits(x, "ltraj"))
        stop("x should be of class \"ltraj\"")
    if (!attr(x, "typeII"))
        stop("x should be of type II (time recorded)")
    res <- lapply(x, function(i) {i$date <- 1:nrow(i);return(i)})
    class(res) <- c("ltraj","list")
    attr(res, "typeII") <- FALSE
    attr(res, "regular") <- FALSE
    res <- rec(res)
    return(res)
}
