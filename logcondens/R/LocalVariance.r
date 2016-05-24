LocalVariance <- function(x, w = NULL, phi){

    n <- length(x)
    if (is.null(w)){p <- rep(1 / n, n)} else {p <-  w / sum(w)}

    m <- sum(p * x)

    dx <- diff(x)
    s2a <- sum(dx * (x[1:(n - 1)] - m) ^ 2 * J10(phi[1:(n - 1)], phi[2:n]))
    s2b <- sum(dx * (x[2:n] - m) ^ 2 * J10(phi[2:n], phi[1:(n - 1)]))
    s2c <- sum(dx ^ 3 * J11(phi[1:(n - 1)], phi[2:n]))

    s2 <- s2a + s2b - s2c

return(s2)
}


