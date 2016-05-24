YP3 <- function (y, d, Z, b1, b2, k, lam, fun) 
{
    yorder <- order(y, -d)
    ysort <- y[yorder]
    dsort <- d[yorder]
    Zsort <- as.matrix(Z[yorder, ])
    gam <- exp(-Zsort %*% cbind(b1, b2))
    n <- length(d)
    S0 <- ((n - 1):0)/n
    fvec <- fun(ysort)
    temp3 <- Haz3(d = dsort, S = S0, gam = gam, lam = lam, fvec = fvec)
    tempT <- Haz3(d = dsort, S = temp3$Su, gam = gam, lam = lam, 
        fvec = fvec)
    for (i in 1:k) {
        tempT <- Haz3(d = dsort, S = tempT$Su, gam = gam, lam = lam, 
            fvec = fvec)
    }
    mu <- sum((tempT$Hazw) * fvec)
    list(d = dsort, Hazw = tempT$Hazw, Survival = tempT$Su, gam = gam, 
        mu = mu)
}
