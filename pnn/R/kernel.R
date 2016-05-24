# D square
# @param Xa A vector describing a new observation
# @param X The training set of observation
ds <- function(Xa, X) {
    value <- (X - Xa) %*% t(X - Xa)
    return(as.numeric(value))
}

# Exponential kernel
# @param Xa A vector describing a new observation
# @param X The training set of observation
# @param sigma The smooth parameter
pattern <- function(Xa, X, sigma) {
    res <- exp( - ds(Xa, X) / (2 * sigma ^ 2) )
    return(as.numeric(res))
}

# Apply kernel over all patterns from category A
# @param Xa A vector describing a new observation
# @param X The training set of observation
# @param sigma The smooth parameter
patterns <- function(Xa, X, sigma)
    apply(Xa, 1, pattern, X, sigma)

# Sum the results of applying the kernel over all patterns
# @param X Pattern from which we have to decide a category. It is a set of measurements represented by a p-dimensional vector
# @param Xa One of the training patterns from category A
# @param sigma Smoothing parameter
fA <- function(Xa, X, sigma) {
    if(missing(Xa)) stop("Xa is missing")
    if(missing(X)) stop("X is missing")
    if(missing(sigma)) stop("sigma is missing")
    p <- length(X) # Dimensionality of measurement space
    m <- length(Xa[,1]) # Total number of training patterns from category A
    f <- 1 /((2 * pi) ^ (p / 2) * sigma ^ p) / m * sum(patterns(Xa, X, sigma)) # Probability density function
    return(f)
}
