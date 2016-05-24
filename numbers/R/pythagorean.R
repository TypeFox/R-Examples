##
##  p y t h a g o r e a n . R  Pythogarean Triples
##


pythagorean_triples <- function(c1, c2) {
    stopifnot(is.numeric(c1),  is.numeric(c2), 
              length(c1) == 1, length(c2) == 1)
    if (!isNatural(c1) || !isNatural(c2) || max(c1, 5) > c2)
        stop("Arguments 'c1', 'c2' must be integers with 'max(c1, 5) <= c2'.")
    c1 <- max(c1, 5)

    M <- c()
    for (cc in c1:c2) {
        m1 <- ceiling(sqrt(cc/2))
        m2 <- floor(sqrt(cc))
        for (m in m1:m2) {
            n <- sqrt(cc - m^2)
            if (floor(n) == ceiling(n)) {
                if (coprime(m, n) && (m-n) %% 2 == 1) {
                    # cat(m, n, "\t", m^2 - n^2, 2*m*n, m^2 + n^2, "\n")
                    M <- rbind(M, c(m, n, m^2 - n^2, 2*m*n, m^2 + n^2))
                }
            }
        }
    }
    return(M)
}
