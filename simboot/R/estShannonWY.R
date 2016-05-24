estShannonWY <- function (x)
{
    Nind <- sum(x)
    estp <- as.numeric(x)/Nind
    estpo0 <- estp[estp > 0]
    Nspecpres <- length(estpo0)
    estH <- (-1) * sum(estpo0 * log(estpo0))
    estHBCo03rd <- estH + (Nspecpres - 1)/(2 * Nind) - (1-sum(1/estpo0))/(12*Nind^2) - sum((1/estpo0)-(1/(estpo0^2)))/(12*Nind^3)
    return(estHBCo03rd)
}
