##
##  c o n t F r a c . R  Continuous Fractions
##


contFrac <- function(x, tol = 1e-6) {
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
