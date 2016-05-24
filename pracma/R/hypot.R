##
##  h y p o t h . R
##


hypot <- function(x, y) {
    if ((length(x) == 0 && is.numeric(y) && length(y) <= 1) ||
        (length(y) == 0 && is.numeric(x) && length(x) <= 1))
            return(c())
    if (!is.numeric(x) && !is.complex(x) || !is.numeric(y) && !is.complex(y))
            stop("Arguments 'x' and 'y' must be numeric or complex.")
    if ((is.vector(x) && is.vector(y) && length(x) != length(y)) ||
        (is.matrix(x) && is.matrix(y) && dim(x) != dim(y)) ||
        (is.vector(x) && is.matrix(y)) || is.matrix(x) && is.vector(y))
            stop("Arguments 'x' and 'y' must be of the same size.")

    x <- abs(x); y <- abs(y)
    m <- pmin(x, y); M <- pmax(x, y)
    ifelse(M == 0, 0, M * sqrt(1 + (m / M)^2))
}
