############################################################################################
defaultmodel <- function (CL, detectfn) {
    if (detectfn %in% c(0:8))
        model <- list(g0 = ~ 1, sigma = ~ 1)
    else if (detectfn %in% c(9))
        model <- list(b0 = ~ 1, b1 = ~ 1)
    else if (detectfn %in% c(10:11))
        model <- list(beta0 = ~ 1, beta1 = ~ 1)
    else ## detectfn %in% c(14:18))
        model <- list(lambda0 = ~ 1, sigma = ~ 1)
    if (!CL) model <- c(list(D = ~1), model)
    model
}
############################################################################################

