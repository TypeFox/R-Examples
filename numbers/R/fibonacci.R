##
##  f i b o n a c c i . R  Fibonacci Sequence
##


fibonacci <- function(n, sequence = FALSE) {
    if (!is.numeric(n) || length(n) != 1 || floor(n) != ceiling(n) || n < 0)
        stop("Argument 'n' must be a single integer >= 0.")
    if (n <= 0) return(0)
    if (n == 1) return(1)
    if (n > 78)
        warning("For 'n > 78' not exactly representable in R as integer.")
    
    if (sequence) {
        if (n == 2) return(c(1, 1))
        fib <- numeric(n)
        fib[1:2] <- c(1, 1)
        for (k in 3:n) {
            fib[k] <- fib[k-1] + fib[k-2]
        }
    } else {
        if (n == 2) return(1)
        f1 <- 1; f2 <- 1
        for (i in 1:(n-2)) {
            t <- f2; f2 <- f1 + f2; f1 <- t
        }
        fib <- f2
    }
    return(fib)
}


catalan <- function(n) {
    if (!is.numeric(n) || length(n) != 1 || floor(n) != ceiling(n) || n < 0)
        stop("Argument 'n' must be a single integer >= 0.")
    if (n >= 30)
        warning("For 'n > 30' this will generate double non-integers.")
    
    if (n <= 1) return(c(1))
    C <- 1
    # this will generate intermediate integers for n <= 30
    for (i in 1:(n-1))
        C <- (n + 1 + i) * C / i

    return(C/n)
}


lucas <- function(n) {
    if (!is.numeric(n) || length(n) != 1 || floor(n) != ceiling(n))
        stop("Argument 'n' must be a single integer >= 0.")
    if (n == 0) return(c(2))
    if (n == 1) return(c(1))
    if (abs(n) > 76)
        warning("For 'n > 76' not exactly representable in R as integer.")

    if (n > 1) {
        l1 <- 2; l2 <- 1
        for (i in 1:(n-1)) {
            t <- l2; l2 <- l1 + l2; l1 <- t
        }
        luc <- l2
    } else {  # n < 0
        luc <- (-1)^n * lucas(-n)
    }
    return(luc)
}


bell <- function(n) {
    stopifnot(is.numeric(n), length(n) == 1)
    if (n < 0 || floor(n) != ceiling(n))
        stop("Argument 'n' must be a whole number greater or equal zero.")
    if (n == 0 || n == 1) return(1)

    B <- Bneu <- numeric(n)
    B[1] <- 1
    for (i in 1:(n-1)) {
        Bneu[1] <- B[i]
        for (j in 2:(i+1)) {
            Bneu[j] <- B[j-1] + Bneu[j-1]
        }
    B <- Bneu
    }
    Bneu[i+1]
}


zeck <- function(n) {
    stopifnot(is.numeric(n))
    if (!isNatural(n) || length(n) != 1)
        stop("Argument 'n' must be an integer.")
    if (n == 1) return(list(fibs = 1, inds = 1))

    Fib <- c(1, 2)
    k <- 2
    f <- 3
    while (f <= n) {
        Fib <- c(Fib, f)
        f <- Fib[k] + f
        k <- k + 1
    }

    fib <- Fib
    K <- c()
    while (n > 0) {
        K <- c(K, k)
        n <- n - fib[k]
        fib <- fib[fib <= n]
        k <- length(fib)
    }

    K <- rev(K)
    return(list(fibs = Fib[K], inds = K+1))
}
