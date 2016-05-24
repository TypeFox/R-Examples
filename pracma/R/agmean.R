##
##  a g m e a n . R  Arithmetic-geometric Mean
##


agmean <- function(a, b) {
    eps <- .Machine$double.eps
    stopifnot(is.numeric(a) || is.complex(a),
              is.numeric(b) || is.complex(b))
    if (is.numeric(a) && any(a < 0) || is.numeric(b) && any(b < 0)) {
        a <- as.complex(a)
        b <- as.complex(b)
    }

    if (length(a) == 1) {
        n <- length(b)
        a <- rep(a, n)
    } else if (length(b) == 1) {
        n <- length(a)
        b <- rep(b, n)
    } else 
        if (length(a) != length(b)) 
            stop("Arguments must have the same length or one has length 1.")

    niter = 0
    while ( any(abs(a-b) >= eps) ) {
        niter = niter + 1
        a1 <- (a + b) / 2
        b1 <- sqrt(a * b)
        if (max(abs(a-a1)) < eps && max(abs(b-b1)) < eps) break
        a <- a1
        b <- b1
    }

    return( list(agm = (a+b)/2, niter = niter, prec = max(abs(b-a))) )
}

