MaximumCost <-
function(x, y, rho.fn, ...) {
    ## Helper function to report the maximum cost produced by the
    ## function rho.fn over the input and output alphabets x and y.
    
    if(class(x) != "matrix") {
        x <- as.matrix(x)
    }
    
    if(class(y) != "matrix") {
        y <- as.matrix(y)
    }
    
    nj <- nrow(x)
    nk <- nrow(y)
    
    max.cost <- -Inf
    
    for(j in 1:nj) {
        xj <- matrix(data = x[j, ], nrow = nk, ncol = ncol(x), byrow = TRUE)
        cost <- rho.fn(xj, y, ...)
        max.cost <- max(max.cost, max(cost))
    }
    for(k in 1:nk) {
        yk <- matrix(data = y[k, ], nrow = nj, ncol = ncol(y), byrow = TRUE)
        cost <- rho.fn(x, yk, ...)
        max.cost <- max(max.cost, max(cost))
    }

    max.cost
}
