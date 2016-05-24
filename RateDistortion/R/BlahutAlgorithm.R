BlahutAlgorithm <-
function(x, px, y, rho.fn, s,
                            eps = 0.001,
                            max.iters = Inf,
                            rho.scale = 1, ...) {
    ## Implementation of the basic Blahut algorithm for solving a
    ## rate-distortion problem at the point s (representing the slope
    ## of the rate-distortion curve).

    if(class(x) != "matrix") {
        x <- as.matrix(x)
    }

    if(class(y) != "matrix") {
        y <- as.matrix(y)
    }

    px <- px / sum(px)
    
    nj <- nrow(x)
    nk <- nrow(y)
    
    q <- 1/nk + vector(mode = "numeric", length = nk)
    
    c <- vector(mode = "numeric", length = nk)
    a <- vector(mode = "numeric", length = nj)
    
    Aj <- vector(mode = "list", length = nj)
    Ak <- vector(mode = "list", length = nk)

    for(j in 1:nj) {
        xj <- matrix(data = x[j, ], nrow = nk, ncol = ncol(x), byrow = TRUE)
        Aj[[j]] <- exp(s * rho.scale * rho.fn(xj, y, ...))
    }
    
    for(k in 1:nk) {
        yk <- matrix(data = y[k, ], nrow = nj, ncol = ncol(y), byrow = TRUE)
        Ak[[k]] <- exp(s * rho.scale * rho.fn(x, yk, ...))
    }

    ##****************************************
    ## Iteration

    termination.reason <- 0
    iters <- 0
    while(TRUE) {
        iters <- iters + 1
        for(j in 1:nj) {
            a[j] <- sum(q * Aj[[j]])
        }

        for(k in 1:nk) {
            c[k] <- sum(px * Ak[[k]] / a)
        }
        
        q <- q * c

        log.c <- log(c)
        TU <- -sum(q * log.c)
        TL <- -max(log.c)

        if((TU - TL) < eps) {
            termination.reason <- 0
            break
        }

        if(iters >= max.iters) {
            termination.reason <- 1
            break
        }
    }

    ##****************************************
    ## Compute R & D

    D <- 0
    for(j in 1:nj) {
        Qj <- Aj[[j]] * q
        Qj <- Qj / sum(Qj)

        xj <- matrix(data = x[j, ], nrow = nk, ncol = ncol(x), byrow = TRUE)
        D <- D + sum(px[j] * Qj * rho.scale * rho.fn(xj, y, ...))
    }

    channel <- list(
        x = x, px = px, y = y,
        rho.fn = rho.fn, s = s,
        eps = eps, iters = iters,
        R = NA, D = D / rho.scale,
        q = q, c = c,
        Aj = Aj, Ak = Ak,
        rho.scale = rho.scale,
        termination = termination.reason)
    class(channel) <- "channel"

    ## Compute the channel information rate
    if(TRUE) {
        ## Method reported by Blahut (1972). Faster than direct computation of mutual information.
        R <- s * D
        for(j in 1:nj) {
            R <- R  - px[j] * log(sum(q * Aj[[j]]))
        }
        R <- R - sum(q * c * log(c))
        R <- R * log2(exp(1)) # Convert from nats to bits
        channel$R <- R
    } else {
        channel$R <- MutualInformation(channel)
    }
    channel
}
