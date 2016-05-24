ChannelDistortion <- function(channel, rho.fn, ...) {
    ## This function returns the expected distortion for a given
    ## channel according to the cost function rho.fn.

    if(class(channel$x) != "matrix") {
        x <- as.matrix(x)
    } else {
        x <- channel$x
    }
    
    if(class(channel$y) != "matrix") {
        y <- as.matrix(y)
    } else {
        y <- channel$y
    }
    
    nj <- nrow(channel$x)
    nk <- nrow(channel$y)
    
    D <- 0
    for(j in 1:nj) {
        py <- ConditionalDistribution(channel, j)
        xj <- matrix(data = x[j, ], nrow = nk, ncol = ncol(x), byrow = TRUE)
        D <- D + channel$px[j] * sum(py$p * rho.fn(xj, y, ...))
    }
    D
}
