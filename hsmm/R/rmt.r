# random values for multivatiate t-distribution
# source: package csampling

rmt <- function (n, df = stop("'df' argument is missing, with no default"), 
    mm = rep(0, mult), cov = diag(rep(1, mult)), mult, is.chol = FALSE) 
{
    if (!missing(cov) && (dim(as.matrix(cov))[1] != dim(as.matrix(cov))[2])) 
        stop("covariance matrix not square")
    if (!missing(mm) && !missing(cov) && (length(mm) != dim(as.matrix(cov))[2])) 
        stop("size mismatch for mean vector and covariance matrix")
    if (missing(mult)) {
        mult <- if (!missing(mm)) 
            length(mm)
        else if (!missing(cov)) 
            dim(as.matrix(cov))[2]
        else 1
    }
    else {
        if (!missing(mm) && (length(mm) != mult)) 
            stop("size mismatch")
        if (!missing(cov) && (dim(as.matrix(cov))[2] != mult)) 
            stop("size mismatch")
    }
    if (mult == 1) 
        return(sqrt(cov) * rt(n, df = df) + mm)
    S <- if (!is.chol) 
        t(chol(cov))
    else cov
    x <- matrix(rnorm(mult * n), nrow = mult, ncol = n, byrow = FALSE)
    y <- rchisq(n, df)
    a <- if (length(n) > 1) 
        S %*% x %*% diag(as.vector(sqrt(df/y))) + matrix(c(rep(mm, 
            n)), nrow = mult, ncol = n, byrow = FALSE)
    else S %*% x * sqrt(df/y) + matrix(c(rep(mm, n)), nrow = mult, 
        ncol = n, byrow = FALSE)
    names(a) <- c()
    a
}
