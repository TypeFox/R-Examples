##
##  p e r m s . R  Permutations
##

perms <- function(a) {
    n <- length(a)
    if (length(a) == 0) return(c())
    else    if (n <= 1) return(matrix(a, 1, 1))

    q <- perms(1:(n-1))  # recursive call
    m <- nrow(q)
    P <- matrix(0, n*m, n)
    P[1:m, ] <- cbind(matrix(n, m, 1), q)

    for (i in (n-1):1) {
        t <- q
        t[t == i] <- n
        P[(m*(n-i)+1):(m*(n-i+1)), ] <- cbind(i*matrix(1, m, 1), t)
    }

    b <- a[c(P)]
    dim(b) <- dim(P)
    return(b)
}


randperm <- function(a, k) {
    n <- length(a)
    if (n == 0 || a[1] == 0) return(c())
    if (n == 1) {
        if (floor(a) != ceiling(a) || a < 1)
            stop("Argument 'a' must be a positive integer.")
        n <- a; a <- 1:a
    }
    if (missing(k)) k <- n
    if (k > n)
        stop("'k' must be smaller or equal to 'a' or length of 'a'.")

    m <- sample(a, size = k, replace = FALSE)
    return(m)
}
