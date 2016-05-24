##
##  h i s t c . R  Histogram Count
##


histc <- function(x, edges) {
    stopifnot(is.numeric(x), is.numeric(edges))

    edges <- c(edges)
    n <- length(edges)
    if (is.unsorted(edges))
        stop("Argument 'edges' must be a monotonically non-decreasing vector.")

    if (length(edges) == 1) {
        bin <- numeric(length(x))
        if (is.matrix(x)) dim(bin) <- c(n, ncol(x))
    	return(list(cnt = 0, bin = bin))
    }

    bin <- numeric(length(x))
    if (is.vector(x)) {
        cnt <- numeric(n)
        for (i in 1:(n-1)) {
            li <- edges[i] <= x & x < edges[i+1]
            cnt[i] <- sum(li)
            bin[li] <- i
        }
        li <- x == edges[n]
        cnt[n] <- sum(li)
        bin[li] <- n

    } else if (is.matrix(x)) {
        cnt <- matrix(0, n, ncol(x))
        for (i in 1:(n-1)) {
            li <- edges[i] <= x & x < edges[i+1]
            cnt[i, ] <- apply(li, 2, sum)
            bin[li] <- i
        }
        li <- x == edges[n]
        cnt[n, ] <- apply(li, 2, sum)
        bin[li] <- n

    } else {
        stop("Argument 'x' must be a numeric vector or matrix.")
    }

    dim(bin) <- dim(x)
    return(list(cnt = cnt, bin = bin))
}


histss <- function(x, n = 100, plotting = FALSE) {
    stopifnot(is.numeric(x), is.numeric(n))
    x <- c(x)
    if (length(n) > 1 || n < 2 || floor(n) != ceiling(n))
        stop("Argument 'n' must be an integer greater than 1.")

    D <- C <- numeric(n-1)
    for (i in 1:(n-1)) {
        D[i] <- diff(range(x)) / (i+1)

        E    <- seq(min(x), max(x), length.out = i+1)
        hp   <- hist(x, breaks = E, plot = FALSE)
        ki   <- hp$counts
        k    <- mean(ki)
        v    <- sum((ki-k)^2) / (i+1)

        C[i] <- (2*k - v) / D[i]^2                  # cost function
    }

    idx  <- which.min(C)
    optD <- D[idx]
    E    <- seq(min(x), max(x), length = idx+1)
    h    <- hist(x, breaks = E, plot = plotting)    # rug(x)

    if (plotting) invisible(h)
    else          return(h)
}
