dg <- function(x, a, b) {
    l <- max(length(x), length (a), length(b))
    a <- rep(a,l)[1:l]
    b <- rep(b,l)[1:l]
    x <- rep(x,l)[1:l]
    dgamma(x, a, b)
}

n <- 10000
x <- seq(0, 10, len = n)
R <- 20
for (i in 1 : R) dg(x, 2, 3)
