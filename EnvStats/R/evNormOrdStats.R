evNormOrdStats <-
function (n = 1, approximate = FALSE) 
{
    if (length(n) > 1 || n != trunc(n) || n < 1) 
        stop("'n' must be a positive integer")
    if (n == 1) 
        return(0)
    if (approximate) {
        Z <- qnorm(ppoints(n, a = 0.375))
    }
    else {
        Z <- numeric(n)
        mid.point <- ifelse(ion <- is.odd(n), (n + 1)/2, n/2)
        if (ion) {
            for (i in 1:(mid.point - 1)) Z[i] <- evNormOrdStatsScalar(i, 
                n)
            Z[(mid.point + 1):n] <- -Z[(mid.point - 1):1]
            Z[mid.point] <- 0
        }
        else {
            for (i in 1:mid.point) Z[i] <- evNormOrdStatsScalar(i, 
                n)
            Z[(mid.point + 1):n] <- -Z[mid.point:1]
        }
    }
    Z
}
