## 
##  s u m a l t . R  Summing alternating series
##


sumalt <- function(f_alt, n) {
    b <- 2^(2*n-1)
    c <- b  # ; s <- 0
    s <- 0.0
    for (k in (n-1):0) {
        t <- f_alt(k)
        s <- s + c*t
        b <- b * (2*k+1) * (k+1) / (2 * (n-k) * (n+k))
        c <- c + b
    }
    s <- s / c
    return(s)
}
