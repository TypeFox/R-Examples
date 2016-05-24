kern <-
function (a, z, lz) 
{
    cat("Calculation of the kernel estimator", fill = TRUE)
    n2 <- length(a)
    band1 <- (sqrt(var(a)) * (n2^(-1/5))) * 1.144
    kernx <- rep(0, lz)
    for (k in 1:lz) {
        zk <- rep(z[k], n2)
        kernxx <- dnorm(((zk - a)/band1), 0, 1)
        kernx[k] <- sum(kernxx)/(n2 * band1)
    }
    return(kernx)
}

