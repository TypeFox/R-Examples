###
### POLYDER.R  Polynom
###

polyder <- function(p, q) {
    if (!missing(q)) {
        if (length(q) == 0) return(0)
        if (!is.numeric(q) && !is.complex(q))
            stop("Arguments must be real or complex vectors or matrices.")
        m <- length(q)
        if (is.matrix(q)) q <- q[1:m]
    } else {
        q <- 1; m <- 1
    }

    if (length(p) == 0) return(0)
    if (!is.numeric(p) && !is.complex(p))
        stop("Argument 'p' must be a real or complex vector or matrix.")
    n <- length(p)
    if (is.matrix(p)) p <- p[1:n]

    # multiply polynomials p an q
    if (n*m <= 1) {
        return(0)
    } else {
        r <- rep(0, n+m-1)
        for (i in seq(along=q)) {
            r <- r + c(rep(0, i-1), p * q[i], rep(0, m-i))
        }
    }

    # case k > 1
    k <- length(r)
    r <- c((k-1):1) * r[1:(k-1)]

    while (r[1] == 0 && length(r) >= 2) {
        r <- r[2:length(r)]
    }
    return(r)
}
