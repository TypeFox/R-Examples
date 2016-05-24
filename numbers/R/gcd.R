##
##  g c d . R  GCD and LCM
##


extGCD <- function(a, b) {
    # The Blankinship method, MAA Mathematical Monthly, Jul-Aug 1963
    stopifnot(is.numeric(a), length(a) == 1, floor(a) == ceiling(a), 
              is.numeric(b), length(b) == 1, floor(b) == ceiling(b))

    sign_ab <- sign(c(a, b))
    A <- matrix(c(abs(c(a, b)), 1, 0, 0, 1), nrow=2, ncol=3)

    while (A[1, 1]*A[2, 1] != 0) {
        if (A[1, 1] > A[2, 1]) {
            m <- A[1, 1] %/% A[2, 1]
            A[1, ] <- A[1, ] - m * A[2, ]
        } else {
            m <- A[2, 1] %/% A[1, 1]
            A[2, ] <- A[2, ] - m * A[1, ]
        }
    }

    if (A[1, 1] == 0)  g <- A[2, ]
    else               g <- A[1, ]

    g[2:3] <- sign_ab * g[2:3]
    return(g)
}


GCD <- function(n, m) {
    stopifnot(is.numeric(n), is.numeric(m))
    if (length(n) != 1 || floor(n) != ceiling(n) ||
        length(m) != 1 || floor(m) != ceiling(m))
        stop("Arguments 'n', 'm' must be integer scalars.")
    if (n == 0 && m == 0) return(0)

    n <- abs(n); m <- abs(m)
    if (m > n) {
        t <- n; n <- m; m <- t
    }
    while (m > 0) {
        t <- n
        n <- m
        m <- t %% m
    }
    return(n)
}


LCM <- function(n, m) {
    stopifnot(is.numeric(n), is.numeric(m))
    if (length(n) != 1 || floor(n) != ceiling(n) ||
        length(m) != 1 || floor(m) != ceiling(m))
        stop("Arguments 'n', 'm' must be integer scalars.")
    if (n == 0 && m == 0) return(0)

    return(n / GCD(n, m) * m)
}


coprime <- function(n, m) {
    stopifnot(is.numeric(n), is.numeric(m))
    if (length(n) != 1 || floor(n) != ceiling(n) ||
        length(m) != 1 || floor(m) != ceiling(m))
        stop("Arguments 'n', 'm' must be integer scalars.")
    if (n == 0 && m == 0) return(FALSE)

    if (GCD(n, m) > 1) FALSE else TRUE
}


mGCD <- function(x) {
    stopifnot(is.numeric(x))
    if (floor(x) != ceiling(x) || length(x) < 2)
        stop("Argument 'x' must be an integer vector of length >= 2.")

    x <- x[x != 0]
    n <- length(x)
    if (n == 0) {
    	g <- 0
    } else if (n == 1) {
    	g <- x
    } else if (n == 2) {
    	g <- GCD(x[1], x[2])
    } else {
        g <- GCD(x[1], x[2])
        for (i in 3:n) {
            g <- GCD(g, x[i])
            if (g == 1) break
        }
    }
    return(g)
}


mLCM <- function(x) {
    stopifnot(is.numeric(x))
    if (floor(x) != ceiling(x) || length(x) < 2)
        stop("Argument 'x' must be an integer vector of length >= 2.")

    x <- x[x != 0]
    n <- length(x)
    if (n == 0) {
    	l <- 0
    } else if (n == 1) {
    	l <- x
    } else if (n == 2) {
    	l <- LCM(x[1], x[2])
    } else {
        l <- LCM(x[1], x[2])
        for (i in 3:n) {
            l <- LCM(l, x[i])
        }
    }
    return(l)
}

