estShannon <- function (x)
{
    Nind <- sum(x)
    estp <- as.numeric(x)/Nind
    estpo0 <- estp[estp > 0]
    pvec <- matrix(estpo0, ncol = 1)
    m1 <- diag(estpo0) - (pvec %*% t(pvec))
    sigma2 <- (t(log(pvec)) %*% m1) %*% log(pvec)
    Nspecpres <- length(estpo0)
    estH <- (-1) * sum(estpo0 * log(estpo0))
    estHBCo03rd <- estH + (Nspecpres - 1)/(2 * Nind) - (1-sum(1/estpo0))/(12*Nind^2) - sum((1/estpo0)-(1/(estpo0^2)))/(12*Nind^3)

    varest <- as.numeric(sigma2)/Nind

    return( list(estimate = estHBCo03rd, varest = varest))
}
