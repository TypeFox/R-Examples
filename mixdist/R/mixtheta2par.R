## last modified June 2002

mixtheta2par <- function(mixtheta, mixpar, constr, mixprop = TRUE) 
{
    k <- nrow(mixpar)
    lentheta <- length(mixtheta)
    if (mixprop) {
        if (constr$conpi == "NONE") {
            lpi <- k - 1
            if (lpi > 0) {
                pi <- mixtheta[1:lpi]
                pi <- c(pi, 1 - sum(pi))
            }
            else pi <- mixpar[, 1]
        }
        else if (constr$conpi == "PFX" & sum(constr$fixpi) < 
            k - 1) {
            pi <- mixpar[, 1]
            lpi <- k - 1 - sum(constr$fixpi)
            lapi <- 1 - sum(mixtheta[1:lpi]) - sum(mixpar[constr$fixpi, 
                1])
            fpi <- c(mixtheta[1:lpi], lapi)
            pi[!constr$fixpi] <- fpi
        }
        else if (constr$conpi == "PFX" & sum(constr$fixpi) >= 
            k - 1) {
            lpi <- 0
            pi <- mixpar[, 1]
        }
    }
    else {
        lpi <- 0
        pi <- mixpar[, 1]
    }
    if (constr$conmu == "NONE") 
        mu <- mixtheta[(lpi + 1):(lpi + k)]
    else if (constr$conmu == "MFX") {
        mu <- mixpar[, 2]
        mu[!constr$fixmu] <- mixtheta[(lpi + 1):(lpi + sum(!constr$fixmu))]
    }
    else if (constr$conmu == "MEQ") 
        mu <- rep(mixtheta[lpi + 1], k)
    else if (constr$conmu == "MES") {
        mu <- mixpar[, 2]
        mu[1:2] <- mixtheta[(lpi + 1):(lpi + 2)]
        if (k >= 3) 
            mu[3:k] <- mu[1] + ((3:k) - 1) * (mu[2] - mu[1])
    }
    else if (constr$conmu == "MGC") {
        mu <- mixpar[, 2]
        mu[1:3] <- mixtheta[(lpi + 1):(lpi + 3)]
        if (k >= 4) 
            mu[4:k] <- mu[1] + ((mu[2] - mu[1])^2) * (1 - ((mu[3] - 
                mu[2])/(mu[2] - mu[1]))^((4:k) - 1))/((mu[2] - 
                mu[1]) - (mu[3] - mu[2]))
    }
    if (constr$consigma == "NONE") 
        sigma <- exp(mixtheta[(lentheta - k + 1):lentheta])
    else if (constr$consigma == "SFX") {
        sigma <- mixpar[, 3]
        sigma[!constr$fixsigma] <- exp(mixtheta[(lentheta - sum(!constr$fixsigma) + 
            1):lentheta])
    }
    else if (constr$consigma == "FCV") 
        sigma <- constr$cov * mu
    else if (constr$consigma == "CCV") {
        cov <- exp(mixtheta[lentheta])/mu[1]
        sigma <- cov * mu
    }
    else if (constr$consigma == "SEQ") 
        sigma <- rep(exp(mixtheta[lentheta]), k)
    else if (constr$consigma == "BINOM") 
        sigma <- sqrt(mu - mu^2/constr$size)
    else if (constr$consigma == "NBINOM") 
        sigma <- sqrt(mu^2/constr$size + mu)
    else if (constr$consigma == "POIS") 
        sigma <- sqrt(mu)
    data.frame(pi = pi, mu = mu, sigma = sigma)
}
