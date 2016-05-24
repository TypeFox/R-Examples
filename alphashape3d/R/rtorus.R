rtorus <-
function (n, r, R, ct = c(0, 0, 0), rotx = NULL) 
{
    if (r < 0 | R < 0 | r > R) {
        stop("The input parameters must satisfy 0 <= r <= R")
    }
    norig <- n
    dif <- R - r
    rsq <- r^2
    sumRr <- r + R
    x <- rep(0, 3)
    while (n > 0) {
        U <- runif(n)
        Rg <- sqrt(U * sumRr^2 + (1 - U) * dif^2)
        Z <- r * (2 * runif(n) - 1)
        phi <- 2 * pi * runif(n)
        aux <- cbind(Rg * cos(phi), Rg * sin(phi), Z)
        pin <- (Rg - R)^2 + Z^2 <= rsq
        if (sum(pin) > 0) {
            x <- rbind(x, aux[pin, ])
        }
        n <- ifelse(sum(pin) < n, n - sum(pin), 0)
    }
    x <- x[-1, ]
    x <- x[1:norig, ]
    x[, 1] <- ct[1] + x[, 1]
    x[, 2] <- ct[2] + x[, 2]
    x[, 3] <- ct[3] + x[, 3]
    if (!is.null(rotx)) {
        aux1 <- cos(rotx) * x[, 2] - sin(rotx) * x[, 3]
        aux2 <- sin(rotx) * x[, 2] + cos(rotx) * x[, 3]
        x[, 2] <- aux1
        x[, 3] <- aux2
    }
    return(x)
}
