NS <- function(param, tm) {
    ## param == c(beta1, beta2, beta3, lambda)
    if (any(tm <= 0))
        stop("all 'tm' must be greater than zero")
    aux <- tm / param[4L]
    param[1L] + param[2L] * ( (1 - exp(-aux)) / aux ) +
        param[3L] * ( (1 - exp(-aux)) / aux - exp(-aux) )
}
NSS <- function(param, tm) {
    ## param = c(beta1, beta2, beta3, beta4, lambda1, lambda2)
    if (any(tm <= 0))
        stop("all 'tm' must be greater than zero")
    gam1 <- tm / param[5L]
    gam2 <- tm / param[6L]
    aux1 <- 1 - exp(-gam1)
    aux2 <- 1 - exp(-gam2)
    param[1L] + param[2L] * (aux1 / gam1) +
        param[3L] * (aux1 / gam1 + aux1 - 1) +
            param[4L] * (aux2 / gam2 + aux2 - 1)
}
