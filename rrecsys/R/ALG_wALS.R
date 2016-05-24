# k: number of features.
# lambda: regularization term.
# scheme: weighting scheme.
# delta: attribute if used in the combined scheme will determine the combination ratio.

# Reference: R. Pan, Y. Zhou, B. Cao, N.  Liu, R. Lukose, M. Scholz, and Q. Yang.  One-Class Collaborative Filtering.


wALS <- function(data, k = 5, lambda = 0.01, scheme = "None!", delta = 0.04) {
    
    
    x <- data@data
    
    colnames(x) <- NULL
    rownames(x) <- NULL
    
    if (ncol(x) < k | nrow(x) < k) 
        stop("Invalid number of features! \n         Less features than the actual number of items or users! Please correct k!")
    
    if (!(scheme %in% c("uni", "uo", "io", "co"))) 
        stop("Specify \"scheme\" attribute. \n    Values: \"uni\" for uniform, \"uo\" for user oriented, \n    \"io\" for item oriented and \"co\" for combinded wighting scheme. ")
    
    # Initializing V with Gaussian random numbers with mean 0 and standard deviation 0.01
    
    V <- rnorm(k * ncol(x), mean = 0, sd = 0.01)
    V <- matrix(V, ncol = k)
    
    
    # Initialize empty U
    U <- matrix(0, nrow = nrow(x), ncol = k)
    
    W <- weightScheme(x, scheme, delta)
    
    ptm <- Sys.time()
    
    
    p <- W * (U %*% t(V))
    resetrrecsysenv()
    
    # updating U and V.
    while (!isConverged(x, p)) {
        # update user features
        U <- lapply(1:nrow(x), function(i) x[i, ] %*% diag(W[i, ]) %*% V %*% ginv(t(V) %*% diag(W[i, ]) %*% V + lambda * sum(W[i, ]) * diag(k)))
        U <- matrix(unlist(U), nrow = nrow(x), byrow = T)
        # update item features
        V <- lapply(1:ncol(x), function(j) x[, j] %*% diag(W[, j]) %*% U %*% ginv(t(U) %*% diag(W[, j]) %*% U + lambda * sum(W[, j]) * diag(k)))
        V <- matrix(unlist(V), nrow = ncol(x), byrow = T)
        
        p <- W * (U %*% t(V))
    }
    
    
    
    cat("Total execution time:", as.numeric(Sys.time() - ptm, units = "secs"), "seconds. \n")
    
    p_wALS <- list(k = k, lambda = lambda, scheme = scheme)
    
    new("wALSclass", alg = "wALS", data = data, factors = list(U = U, V = V), weightScheme = W, parameters = p_wALS)
}

p_wALS <- list(k = 10, lambda = 0.01, scheme = "None!", delta = 0.04)
rrecsysRegistry$set_entry(alg = "wALS", 
                          fun = wALS, 
                          description = "Weighted Alternating Least Squares.", 
                          reference = "R. Pan, Y. Zhou, B. Cao, N.  Liu, R. Lukose, M. Scholz, and Q. Yang.  One-Class Collaborative Filtering.",
                          parameters = p_wALS) 
