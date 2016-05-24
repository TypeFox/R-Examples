NSf <- function(lambda, tm) {
    aux <- tm/lambda
    X <- array(1, dim = c(length(tm),3L))
    X[ ,2L] <- (1 - exp(-aux))/aux
    X[ ,3L] <- ((1 - exp(-aux))/aux) - exp(-aux)
    X
}

NSSf <-function(lambda1, lambda2, tm) {
    aux1 <- tm/lambda1
    aux2 <- tm/lambda2
    X <- array(1, dim = c(length(tm),4L))
    X[ ,2L] <- (1 - exp(-aux1))/aux1
    X[ ,3L] <- (1 - exp(-aux1))/aux1 - exp(-aux1)
    X[ ,4L] <- (1 - exp(-aux2))/aux2 - exp(-aux2)
    X
}
