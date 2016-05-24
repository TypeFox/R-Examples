##
##  k n a p s a c k . R  Knapsack Problems
##


knapsack <- function(w, p, cap) {
    n <- length(w)

    # initiate arrays
    x <- logical(n)
    F <- matrix(0, nrow = cap + 1, ncol = n)
    G <- matrix(0, nrow = cap + 1, ncol = 1)

    # forwarding part
    for (k in 1:n) {
        F[, k] <- G
        H <- c(numeric(w[k]), G[1:(cap + 1 - w[k]), 1] + p[k])
        G <- pmax(G, H)
    }
    fmax <- G[cap + 1, 1]

    # backtracking part
    f <- fmax
    j <- cap + 1
    for (k in n:1) {
        if (F[j, k] < f) {
            x[k] <- TRUE
            j <- j - w[k]
            f <- F[j, k]
        }
    }

    inds <- which(x)
    wght <- sum(w[inds])
    prof <- sum(p[inds])
    return(list(capacity = wght, profit = prof, indices = inds))
}
