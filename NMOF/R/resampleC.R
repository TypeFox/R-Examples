resampleC <- function(..., size, cormat) {
    n <- length(data <- list(...))
    m <- sapply(data, NROW)
    ns <- seq_len(n)
    defnames <- paste("var", ns, sep = "")
    nm <- names(data)
    if (is.null(nm))
        nm <- defnames
    else if (any(hasname <- nzchar(nm)))
        defnames[hasname] <- nm[hasname]
    data <- lapply(data, sort)
    X <- array(rnorm(size*n), dim = c(size, n))
    cormat <- 2*sin(cormat*pi/6)
    ee <- eigen(cormat)
    if (any(ee$values < 0))
        warning("'cormat' has eigenvalues smaller than zero")
    U <- pnorm(X %*% t(ee$vectors %*% diag(sqrt(ee$values))))
    ii <- ceiling(U %*% diag(m))

    results <- array(NA, dim = c(size,n))
    colnames(results) <- defnames
    for (i in seq_len(n))
        results[ ,i] <- data[[i]][ii[ ,i]]
    results
}
