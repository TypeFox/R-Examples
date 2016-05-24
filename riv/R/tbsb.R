Tbsb <- function(c, p) {
    ksiint <- function(s) {
        (2^s) * gamma(s + p/2) * pgamma(c^2/2, s + p/2)/gamma(p/2)
    }
    
    y1 <- ksiint(1) * 3/c - ksiint(2) * 3/(c^3) + ksiint(3)/(c^5)
    y2 <- c * (1 - pchisq(c^2, p))

    y1 + y2
}
