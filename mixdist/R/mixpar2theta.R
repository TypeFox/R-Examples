## last modified June 2002

mixpar2theta <- function(mixpar, constr, mixprop = TRUE) 
{
    k <- nrow(mixpar)
    if (constr$conmu == "NONE") 
        mu.e <- mixpar[, 2]
    else if (constr$conmu == "MFX") 
        mu.e <- mixpar[!constr$fixmu, 2]
    else if (constr$conmu == "MEQ") 
        mu.e <- mixpar[1, 2]
    else if (constr$conmu == "MES") 
        mu.e <- mixpar[1:2, 2]
    else if (constr$conmu == "MGC") 
        mu.e <- mixpar[1:3, 2]
    if (constr$consigma == "NONE") 
        sigma.e <- log(mixpar[, 3])
    else if (constr$consigma == "SFX") 
        sigma.e <- log(mixpar[!constr$fixsigma, 3])
    else if (!is.na(match(constr$consigma, c("FCV", "BINOM", 
        "NBINOM", "POIS")))) 
        sigma.e <- NULL
    else if (!is.na(match(constr$consigma, c("CCV", "SEQ")))) 
        sigma.e <- log(mixpar[1, 3])
    if (mixprop) {
        if (constr$conpi == "NONE") 
            pi.e <- mixpar[-k, 1]
        else if (constr$conpi == "PFX" & sum(constr$fixpi) < 
            k - 1) {
            pi.e <- mixpar[!constr$fixpi, 1]
            lpi <- length(pi.e)
            pi.e <- pi.e[-lpi]
        }
        else if (constr$conpi == "PFX" & sum(constr$fixpi) >= 
            k - 1) 
            pi.e <- NULL
    }
    else pi.e <- NULL
    c(pi.e, mu.e, sigma.e)
}
