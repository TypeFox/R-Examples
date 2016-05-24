FPRP <- function (a, b, pi0, ORlist, logscale = FALSE)
{
    if (logscale) {
        thetahat <- a
        V <- b
    }
    else {
        thetahat <- log(a)
        V <- ((log(b) - log(a))/1.96)^2
    }
    s <- sqrt(V)
    if (thetahat > 0) {
        p <- 2 * (1 - pnorm(thetahat/s))
        q <- qnorm(1 - p/2)
        b1 <- 1 - pnorm(q - log(ORlist)/s)
    }
    else {
        p <- 2 * pnorm(thetahat/s)
        q <- qnorm(1 - p/2)
        b1 <- pnorm(-q - log(ORlist)/s)
    }
    FPRP <- t(p * (1 - pi0)/(p * (1 - pi0) + b1 %o% pi0))
    row.names(FPRP) <- pi0
    colnames(FPRP) <- ORlist
    FNRP <- 1/(1 + (1 - pi0)/pi0 %o% ((1 - p)/(1 - b1)))
    row.names(FNRP) <- pi0
    colnames(FNRP) <- ORlist
    invisible(list(p = p, power = b1, FPRP = FPRP, FNRP = FNRP))
}

# 2-8-2007
