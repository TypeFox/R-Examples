gauher <-
function (n) {
    m <- trunc((n + 1)/2)
    x <- w <- rep(-1, n)
    for (i in seq_len(m)) {
        z <- if (i == 1) {
            sqrt(2*n + 1) - 1.85575 * (2*n + 1)^(-0.16667)
        } else if (i == 2) {
            z - 1.14 * n^0.426 / z
        } else if (i == 3) {
            1.86 * z - 0.86 * x[1]
        } else if (i == 4) {
            1.91 * z - 0.91 * x[2]
        } else {
            2*z - x[i - 2]
        }
        for (its in seq_len(10)) {
            p1 <- 0.751125544464943
            p2 <- 0
            for (j in seq_len(n)) {
                p3 <- p2
                p2 <- p1
                p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
            }
            pp <- sqrt(2*n) * p2
            z1 <- z
            z <- z1 - p1/pp
            if (abs(z - z1) <= 3e-14) 
                break
        }
        x[i] <- z
        x[n + 1 - i] <- -z
        w[i] <- 2 / (pp * pp)
        w[n + 1 - i] <- w[i]
    }
    list(x = x, w = w)
}
