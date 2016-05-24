rmnorm <- function (n = 1, mean = rep(0, d), varcov) 
{
    varcov <- as.matrix(varcov)
    d <- ncol(varcov)
    z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
    y <- t(mean + t(z))
    y
}