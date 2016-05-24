##
##  a c c u m a r r a y . R  Accumulate Vector Elements
##


uniq <- function(a, first = FALSE) {
    if (length(a) == 0)
        return(list(b = c(), m = c(), n = c()))
    if (!is.numeric(a) || !is.vector(a))
        stop("Argument 'a' must be a numeric vector.")

    la <- length(a); n <- numeric(la)
    u  <- unique(a)
    lu <- length(u); m <- numeric(lu)

    mima <- if (first) min else max

    for (i in 1:lu) {
        w <- which(a == u[i])
        m[i] <- mima(w)
        n[w] <- i 
    }

    return(list(b = u, m = m, n = n))
}


accumarray <- function(subs, val, sz = NULL, func = sum, fillval = 0) {
    stopifnot(is.numeric(subs), is.numeric(val))
    subs <- floor(subs)
    val <- c(val)
    if (any(subs < 1))
        stop("Argument 'subs' must be a matrix of integer indices.")

    matrix_p <- TRUE
    if (is.vector(subs)) {
        subs <- as.matrix(subs)
        matrix_p <- FALSE
    }

    n <- nrow(subs); m <- ncol(subs)
    if (length(val) < n)
        stop("Length of 'vals' must not be smaller than no. of rows of 'subs'.")

    dm <- apply(subs, 2, max)
    if (!is.null(sz)) {
        if (length(sz) != ncol(subs) || any(sz < dm))
            stop("Argument 'sz' does not fit with 'subs'.")
        dm <- sz
    }

    if (m == 1) {
        A <- rep(fillval, dm)
        for (i in unique(subs)) {
            A[i] <- func(val[subs == i])
        }
        if (matrix_p) A <- as.matrix(A)

    } else {
        cm <- cumprod(dm[1:(m-1)])
        A <- array(fillval, dim = dm)

        K <- numeric(n)
        for (i in 1:n) {
            K[i] <- subs[i, 1] +  sum(cm * (subs[i, 2:m]-1))
        }
        for (i in unique(K)) {
            A[i] <- func(val[K == i])
        }
    }
    return(A)
}
