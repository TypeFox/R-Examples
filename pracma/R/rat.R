##
##  r a t . R  Continuous Fractions
##


.contfrac <- function(x, tol = 1e-6) {
    if (!is.numeric(x) || is.matrix(x))
        stop("Argument 'x' must be a numeric scalar or vector.")

    if (length(x) > 1) {
        # Compute value of a continuous fraction
        n <- length(x)
        B <- diag(1, 2)
        for (i in seq(along=x)) {
            B <- B %*% matrix(c(x[i], 1, 1, 0), 2, 2)
        }
        return(B[1,1]/B[2,1])

    } else {
        # Generate the continuous fraction of a value
        sgnx <- sign(x)
        x <- abs(x)

        b <- floor(x)
        k <- b
        r <- x - b
        B <- matrix(c(b, 1, 1, 0), 2, 2)
        while ( abs(x - B[1,1]/B[2,1])  > tol) {
            b <- floor(1/r)
            k <- c(k, b)
            r <- 1/r - b
            B <- B %*% matrix(c(b, 1, 1, 0), 2, 2)
        }
        return(list(cf = sgnx * k, rat = c(sgnx*B[1,1], B[2,1]),
                    prec = abs(x - B[1,1]/B[2,1])))
    }
}


rat <- function(x, tol = 1e-6) {
    if (length(x) == 0)
        return(c())
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric vector.")
    xs <- c(x)

    n <- length(xs)
    R <- character(n)
    for (i in 1:n) {
        x <- xs[i]
        k <- .contfrac(x, tol = tol)$cf
        if (length(k) >= 1) {
            cf <- paste("[ ", k[1], sep="")
        }
        if (length(k) >= 2) {
            cf <- paste(cf, "; ", k[2], sep="")
        }
        if (length(k) >= 3) {
            for (j in 3:length(k)) {
                cf <- paste(cf, ", ", k[j], sep="")
            }
            cf <- paste(cf, "]", sep="")
        }
        R[i] <- cf
    }
    return(R)
}

rats <- function(x, tol = 1e-6) {
    if (length(x) == 0)
        return(c())
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric vector.")
    xs <- c(x)

    n <- length(xs)
    R <- numeric(n)
    for (i in 1:n) {
        x <- xs[i]
        k <- .contfrac(x, tol = tol)$rat
        cf <- paste(k[1], "/", k[2], sep="")
        cat(cf, "\n")
        R[i] <- k[1]/k[2]
    }
    invisible(R)    
}
