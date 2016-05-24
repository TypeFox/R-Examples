##### internal functions from Dr. Tibshirani's software package GSA

cor.func <-
function (x, y ) {
    n <- length(y)
    xbar <- x %*% rep(1/n, n)
    sxx <- ((x - as.vector(xbar))^2) %*% rep(1, n)
    sxy <- (x - as.vector(xbar)) %*% (y - mean(y))
    syy <- sum((y - mean(y))^2)
    numer <- sxy/sxx
    sd <- sqrt((syy/sxx - numer^2)/(n - 2))

    tt <- numer/sd

    return(list(tt = tt, numer = numer, sd = sd))
}
