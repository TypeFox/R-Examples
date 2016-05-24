##
##  c o m b s . R  Combinations
##


combs <- function(a, m){
    n <- length(a)
    if (length(a) == 0 || m <= 0) return(c())
    if (m >= n) return(a)
    if (m <= 1) return(matrix(a, n, 1))

    v <- c(a)
    P <- c()
    for (k in 1:(n-m)) {
        Q <- combs(v[(k+1):n], m-1)
        P <- rbind(P, cbind(v[k], Q))
    }
    k <- n-m+1
    Q <- combs(v[(k+1):n], m-1)
    P <- rbind(P, c(v[k], Q))

    b <- a[c(P)]
    dim(b) <- dim(P)
    return(P)
}

randcomb <- function(a, m) {
    n <- length(a)
    if (n == 0) return(c())


    m <- sample(1:n, size = m, replace = FALSE)
    return(a[m])
}
