ds <- function(Xa, X) {
    value <- (X - Xa) %*% t(X - Xa)
    return(as.numeric(value))
}

pattern <- function(Xa, X, sigma) {
    res <- exp( - ds(Xa, X) / (2 * sigma ^ 2) )
    return(as.numeric(res))
}

patterns <- function(Xa, X, sigma)
    apply(Xa, 1, pattern, X, sigma)

Y <- function(Xa, Ya, X, sigma) {
    p <- length(X) # Dimensionality of measurement space
    m <- length(Xa[,1]) # Total number of training patterns from category A
    patterns1 <- patterns(Xa, X, sigma)
    f <- sum(Ya * patterns1) / sum(patterns1)
    return(f)
}
