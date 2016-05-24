# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


rmvtnorm = function (n = 1, k = 2, means = rep(0, k), sigma = diag(k)){
    if (length(sigma) == 1) if (sigma >= -1 & sigma <= 1) 
            sigma = matrix(c(1, sigma, sigma, 1), 2, 2)
    else stop("Correlation magnitude must be less than or equal to 1.")
    
    if (!identical(sigma, t(sigma))) 
        stop("Inappropriate covariance matrix specified.")
    
    if (nrow(sigma) != k | ncol(sigma) != k | length(means) != k) 
        stop("Incorrect covariance matrix dimensions.")
    
    eig = eigen(sigma, symmetric = TRUE)
    if (sum(eig$values < 0) > 0) 
        stop("Inappropriate covariance matrix specified.")
    
    A = eig$vectors %*% diag(sqrt(eig$values))
    x = matrix(rnorm(n * length(means), 0, 1), k, n)
    output = t(means + A %*% x)
    return(output)
}
