CombineChisquare <- function (Xmat, alpha = 0.05, varest = "oneovern")
{
    m <- nrow(Xmat)
    fvec <- numeric(m)
    varvec <- numeric(m)
    for (i in 1:nrow(Xmat)) {
        fhat <- HWf(Xmat[i, ])
        pa <- af(Xmat[i, ])
        n <- sum(Xmat[i, ])
        fvec[i] <- fhat
        varf <- switch(varest, oneovern = 1/n,
                               bailey  = ((1-fhat)^2)*(1-2*fhat)/n + fhat*(1-fhat)*(2-fhat)/(2*n*pa*(1-pa)),
                               stop("unknown option for parameter varest"))
        varvec[i] <- varf
    }
    MeanTheta <- mean(fvec)
    W <- mean(varvec)
    B <- var(fvec)
    T <- W + (m + 1) * B/m
    gamma <- (1 + 1/m) * B/T
    v <- (m - 1) * (1 + m * W/((m + 1) * B))^2
#    stat <- abs(MeanTheta)/sqrt(T)
    stat <- MeanTheta/sqrt(T)
    pvalimp <- 2 * pt(abs(stat), v, lower.tail = FALSE)
    pvalimpd <- pt(stat, v, lower.tail = FALSE) # H1: f > 0
    pvalimpe <- pt(stat, v, lower.tail = TRUE)  # H1: f < 0
    r <- (1 + 1/m) * B/W
    lambda <- (r + 2/(v + 3))/(r + 1)
    gamma <- (1 + 1/m) * B/T
    gammaalt <- (r + 2/(v + 3))/(r + 1)
    fhatimp <- MeanTheta
    llf <- fhatimp - qt(1 - alpha/2, v) * sqrt(T)
    ulf <- fhatimp + qt(1 - alpha/2, v) * sqrt(T)
    return(list(fhatimp = fhatimp, pvalimp = pvalimp, pvalimpd = pvalimpd, pvalimpe = pvalimpe, r = r,
        lambda = lambda, llf = llf, ulf = ulf, fvec = fvec, varvec = varvec,
        gammaalt = gammaalt))
}
