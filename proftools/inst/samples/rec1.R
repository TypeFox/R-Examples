f <- function(n) {
    n <- n - 1
    g(n)
}

g <- function(n) {
    if (n > 0)
        f(n)
    else
        0
}

h <- function(n) {
    for (i in 1 : n)
        f(100)
}

k <- function(n)
    for (i in 1:n) {
        f(100)
        g(100)
    }

system.time(k(25000))
