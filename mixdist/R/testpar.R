## last modified June 2002

testpar <- function(mixpar, dist, constr) 
{
    parvalid <- TRUE
    k <- nrow(mixpar)
    pi <- mixpar[, 1]
    mu <- mixpar[, 2]
    sigma <- mixpar[, 3]
    if (sum(is.na(pi)) > 0 | sum(is.na(mu)) > 0 | sum(is.na(sigma)) > 
        0) 
        parvalid <- FALSE
    else if (sum(is.nan(pi)) > 0 | sum(is.nan(mu)) > 0 | sum(is.nan(sigma)) > 
        0) 
        parvalid <- FALSE
    else if (sum(is.infinite(pi)) > 0 | sum(is.infinite(mu)) > 
        0 | sum(is.infinite(sigma)) > 0) 
        parvalid <- FALSE
    else if (sum(pi < 0) > 0 | sum(pi > 1) > 0) 
        parvalid <- FALSE
    else if (sum(sigma <= 0) > 0) 
        parvalid <- FALSE
    else if (mu[1] <= 0 & !is.na(match(constr$consigma, c("FCV", 
        "CCV", "BINOM", "NBINOM", "POIS")))) 
        parvalid <- FALSE
    else if (constr$consigma == "MGC" & mu[3] - mu[2] > mu[2] - 
        mu[1]) 
        parvalid <- FALSE
    else if (mu[1] <= 0 & !is.na(match(dist, c("lnorm", "gamma", 
        "weibull")))) 
        parvalid <- FALSE
    else if (mu[1] < 0 & !is.na(match(dist, c("binom", "nbinom", 
        "pois")))) 
        parvalid <- FALSE
    else if (dist == "binom" & constr$consigma != "BINOM" & sum(mu <= 
        sigma^2) > 0) 
        parvalid <- FALSE
    else if (dist == "nbinom" & constr$consigma != "NBINOM" & 
        sum(sigma^2 <= mu) > 0) 
        parvalid <- FALSE
    else {
        if (k > 1) {
            if (sum(mu[-k] - mu[-1] > 0) > 0) 
                parvalid <- FALSE
            if (sum(mu[-k] - mu[-1] == 0 & sigma[-k] - sigma[-1] == 
                0) > 0) 
                parvalid <- FALSE
        }
    }
    parvalid
}
