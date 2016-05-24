"ARToPacf" <- function(phi) {
    phik <- phi
    L <- length(phi)
    if (L == 0) 
        return(numeric(0))
    pi <- numeric(L)
    for (k in 1:L) {
        LL <- L + 1 - k
        pi[L + 1 - k] <- a <- phik[LL]
        phikp1 <- phik[-LL]
        if (abs(a) == 1) 
            break
        phik <- (phikp1 + a * rev(phikp1))/(1 - a^2)
    }
    pi
}

"PacfToAR" <- function(pi) {
    L <- length(pi)
    if (L == 0) 
        return(numeric(0))
    if (L == 1) 
        return(pi)
    phik <- pi[1]
    for (k in 2:L) {
        phikm1 <- phik
        phik <- c(phikm1 - pi[k] * rev(phikm1), pi[k])
    }
    phik
} 
