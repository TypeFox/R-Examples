# NSf.R -- version 2010-12-29
NSf <- function(lambda,tm) {
    aux <- tm / lambda
    ll <- length(tm)
    Y <- array(1, dim = c(ll,3L))
    Y[ ,2L] <- ((1 - exp(-aux)) / aux)
    Y[ ,3L] <- (((1 - exp(-aux)) / aux) - exp(-aux))
    Y
}
NSSf <- function(lambda1,lambda2,tm) {
    aux1 <- tm / lambda1
    aux2 <- tm / lambda2
    ll <- length(tm)
    Y <- array(1, dim = c(ll,4L))
    Y[ ,2L] <- ((1 - exp(-aux1)) / aux1)
    Y[ ,3L] <- (((1 - exp(-aux1)) / aux1) - exp(-aux1))
    Y[ ,4L] <- (((1 - exp(-aux2)) / aux2) - exp(-aux2))
    Y
}

