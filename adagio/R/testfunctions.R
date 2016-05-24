##
##  t e s t f u n c t i o n s . R
##


#-- Rosenbrock's non-convex performance test function ------------------------
fnRosenbrock <- function(x) {
    n <- length(x)
    x1 <- x[2:n]
    x2 <- x[1:(n-1)]
    sum(100*(x1 - x2^2)^2 + (1 - x2)^2)
}

grRosenbrock <- function(x) {
    n <- length(x)
    g <- rep(NA, n)
    g[1] <- 2*(x[1] - 1) + 400*x[1] * (x[1]^2 - x[2])
    if (n > 2) {
        ii <- 2:(n-1)
        g[ii] <- 2*(x[ii]-1) + 400*x[ii]*(x[ii]^2-x[ii+1]) + 200*(x[ii]-x[ii-1]^2) 
    }
    g[n] <- 200*(x[n] - x[n-1]^2)
    g
}


#-- Rastrigin's function -----------------------------------------------------
fnRastrigin <- function(x) {
    n <- length(x)
    10*n + sum(x^2 - 10*cos(2*pi*x))
}

grRastrigin <- function(x) {
    2 * x + 20 * pi * sin(2 * pi * x)
}


#-- Nesterov's smooth Chebyshev-Rosenbrock function --------------------------
fnNesterov <- function(x) {
    n <- length(x)
    f <- (1 - x[1])^2 / 4
    for (i in 1:(n-1)) {
        f <- f + (1 + x[i+1] - 2*x[i]^2)^2
    }
    f
}

grNesterov <- function(x) {
    n <- length(x)
    g <- rep(NA, n)
    g[1] <- (x[1] - 1) / 2
    for (i in 1:(n-1)) {
        r = 1 + x[i+1] - 2*x[i]^2
        g[i+1] <- g[i+1] + 2*r
        g[i] <- g[i] - 8*x[i]*r
    }
    g
}


#-- Nesterov's non-smooth Chebyshev-Rosenbrock functions ---------------------
fnNesterov1 <- function(x) {
    n <- length(x)
    f <- (1 - x[1])^2 / 4
    for (i in 1:(n-1)) {
        f <- f + abs(1 + x[i+1] - 2*x[i]^2)
    }
    f
}


#-- Non-smooth test function of Hald and Madsen ------------------------------
# xmin = (0.99987763,  0.25358844, -0.74660757,  0.24520150, -0.03749029 )
# fmin = 0.0001223713
initHald <- function(x) {
    stopifnot(is.numeric(x), length(x) == 5)
    i <- 1:21
    t <- -1 + (i - 1)/10
    (x[1] + x[2] * t) / ( 1 + x[3]*t + x[4]*t^2 + x[5]*t^3 ) - exp(t)
}

fnHald <- function(x) {
    f <- initHald(x)
    max(abs(f))
}

grHald <- function(x) {
    g <- rep(NA, 5)
    f <- initHald(x)
    k <- which.max(abs(f))
    s1 <- sign(f[k])
    t <- -1.0 + (k-1) * 0.1
    a <-  1 + x[3]*t + x[4]*t^2 + x[5]*t^3
    b <- x[1] + x[2] * t
    s2 <- sign(t * a * b)
    if (a == 0) return(g)
    g[1:2] <- s1 * c(1, t) / a
    g[3:5] <- s2 * c(t, t^2, t^3) * (x[1] + x[2] * t) / a^2
    g
}


#-- Shor's piecewise quadratic function --------------------------------------
# starting value c(1,1,1,1,1) in [0, 2]^5
# xmin = (1.1243510, 0.9794616, 1.4777077, 0.9202335, 1.1242916)
# fmin = 22.6001621
fnShor <- function(x) {
    stopifnot(is.numeric(x), length(x) == 5)
    x <- as.matrix(c(x))
    A <- matrix(
    c(0,  2,  1,  1,  3,  0,  1,  1,  0,  1,
      0,  1,  2,  4,  2,  2,  1,  0,  0,  1,
      0,  1,  1,  1,  1,  1,  1,  1,  2,  2,
      0,  1,  1,  2,  0,  0,  1,  2,  1,  0,
      0,  3,  2,  2,  1,  1,  1,  1,  0,  0),
    5, 10, byrow = TRUE)
    b <- as.matrix(c(1, 5, 10, 2, 4, 3, 1.7, 2.5, 6, 4.5))
    f <- 0
    for (i in 1:10) {
        d <- b[i] * sum((x - A[, i])^2)         # d[i] <- ...
        if (d > f) f <- d                       # f <- max(d)
    }
    f
}

grShor <- function(x) {
    stopifnot(is.numeric(x), length(x) == 5)
    x <- as.matrix(c(x))
    A <- matrix(
    c(0,  2,  1,  1,  3,  0,  1,  1,  0,  1,
      0,  1,  2,  4,  2,  2,  1,  0,  0,  1,
      0,  1,  1,  1,  1,  1,  1,  1,  2,  2,
      0,  1,  1,  2,  0,  0,  1,  2,  1,  0,
      0,  3,  2,  2,  1,  1,  1,  1,  0,  0),
    5, 10, byrow = TRUE)
    b <- as.matrix(c(1, 5, 10, 2, 4, 3, 1.7, 2.5, 6, 4.5))
    f <- 0
    for (i in 1:10) {
        d <- b[i] * sum((x - A[, i])^2)
        if (d > f) {
            f <- d
            k <- i
        }
    }
    c(2 * b[k] * (x - A[, k]))
}
