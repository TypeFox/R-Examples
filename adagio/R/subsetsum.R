##
##  s u b s e t s u m . R  Subset Sum Algorithm
##


# Assume S decreasing, no elements > t, total sum >= t
subsetsum <- function(S, t, method = "greedy") {
    stopifnot(is.numeric(S), is.numeric(t))
    if (length(t) > 1) {
        t <- t[1]
        warning("Length of 't' must be 1; will take only first component.")
    }
    if (any(floor(S) != ceiling(S)) || any(S <= 0) ||
        floor(t) != ceiling(t)      || t <= 0) 
        stop("Arguments 'S' and 't' must be positive integer vectors.")

    if (any(S >= t))
        stop("No element of 'S' shall be greater or equal to 't'.")
    if (sum(S) < t) {
        warning("Total sum of 'S' is smaller than 't'; no solution possible.")
        return(NA)
    }

    method <- pmatch(method, c("greedy", "dynamic"))
    if (is.na(method))
        stop("The 'method' must be one of 'greedy' or 'dynamic'.")
    method <- c("greedy", "dynamic")[method]

    n <- length(S)
    inds <- NULL

    if (method == "greedy") {
        L <- c(0)
        for (i in 1:n) {
            L <- unique(c(L, L+S[i]))
            L <- L[L <= t]
            if (max(L) == t) {
                inds <- c(i)
                t <- t - S[i]
                while (t > 0) {
                    K <- c(0)
                    for (j in 1:n) {
                        K <- unique(c(K, K+S[j]))
                        K <- K[K <= t]
                        if (max(K) == t) break
                    }
                    inds <- c(inds, j)
                    t <- t - S[j]
                }
                break
            }
        }
        if (length(inds) != 0) inds <- sort(inds)

    } else {  # if (method == "dynamic")
        x <- logical(n)
        F <- numeric(t + 1)
        G <- logical(t + 1)
        G[1] <- TRUE
        for (k in 1:n) {
            H <- c(logical(S[k]), G[1:(t + 1 - S[k])])
            H <- (G < H)
            j <- which(H)
            F[j] <- k
            G[j] <- TRUE
            if (G[t + 1]) break
        }
        wch <- which(G)
        j <- wch[length(wch)]
        fmax <- j - 1
        while (j > 1) {
            k <- F[j]
            x[k] <- TRUE
            j <- j - S[k]
        }
        inds <- which(x)
    }

    return(list(val = sum(S[inds]), inds = inds))
}
