slices <- function(from, to, n) {
    if (!is.numeric(from) && length(from) != 1) {
        stop("'from' must be a numeric value")
    }
    if (!is.numeric(to) && length(to) != 1) {
        stop("'to' must be a numeric value")
    }
    if (!is.numeric(n) && length(n) != 1) {
        stop("'n' must be a numeric value")
    }
    if (n < 2) {
        stop("'n' must be at least 2")
    }
    res <- 0:(n-1) / (n-1) * (to - from) + from
    return(res)
}
