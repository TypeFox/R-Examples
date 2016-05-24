##
##  r a b i n . R  Rabin-Miller primality test
##


miller_rabin <- function(n) {
    stopifnot(is.numeric(n), length(n) == 1)
    if (floor(n) != ceiling(n) || n <= 0) {
        stop("Argument 'n' must be a natural number.")
    }
    if (n < 2) { return(FALSE)
    } else if (n == 2) {return(TRUE)
    } else if (n %% 2 == 0) {return(FALSE)
    } else if (n >= 2^53) {
        stop("Argument 'n' too large to be handled as integer.")
    }

    if (!requireNamespace("gmp", quietly = TRUE)) {
        stop("Package 'gmp' needed: Please install separately.", call. = FALSE)
    }


    # define witnesses (see Miller-Rabin test at Wolfram MathWorld)
    if (n < 2047) {
        N <- c(2)
    } else if (n < 1373653) {
        N <- c(2, 3)
    } else if (n < 25326001) {
        N <- c(2, 3, 5)
    } else if (n < 3215031751) {
        N <- c(2, 3, 5, 7)
    } else if (n < 2152302898747) {
        N <- c(2, 3, 5, 7, 11)
    } else if (n < 3474749660383) {
        N <- c(2, 3, 5, 7, 11, 13)
    } else if (n < 341550071728321) {
        N <- c(2, 3, 5, 7, 11, 13, 17)
    } else {
        # warning("For 'n > 341550071728321' this test is only probabilistic.")
        N <- c(2, 3, 5, 7, 11, 13, 17, 23)
    }

    # extract powers of 2
    m <- n-1; s <- 0
    while ((m %% 2) == 0) {
        s <- s + 1
        m <- m / 2
    }

    # check witnesses
    for (a in N) {
        if (a >= n) break
        # x <- 1
        # for (i in 1:m) {
        #     if(a*x >= 2^53) stop("Error: Integer arithmetic")
        #     x <- (a*x) %% n
        # }
        x <- gmp::powm(a, m, n)
        if (x == 1) next
        t <- s
        while (x != (n-1)) {
            t <- t-1
            if (t <= 0) return(FALSE)
            x = ((x*x) %% n)
            if (x == 1) return(FALSE)
        }
    }
    # Lucas-Lehmer test
    return(TRUE)
}
