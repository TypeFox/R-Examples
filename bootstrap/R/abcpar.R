"abcpar" <- function(y, tt, S, etahat, mu, n=rep(1,length(y)), lambda =
                     0.001, alpha = c(0.025, 0.05, 0.1, 0.16))
{
    call <- match.call()
    syscall <- sys.call()
    p <- length(y)
    I <- diag(p)
    ## calculate thetahat,ehat,dhat, and sighat
    thetahat <- tt(y)
    ehat <- numeric()
    for(j in 1:p) {
        lam <- lambda * S[j, j]^0.5
        delta <- I[, j]
        ehat[j] <- (tt(y + lam * delta) - tt(y - lam * delta))/(2 * lam)
    }
    dhat <- as.vector(S %*% ehat)
    sighat <- sqrt(ehat %*% S %*% ehat)
    ##  calculate acceleration a
    lam <- lambda/sighat
    a0 <- sum(ehat * mu(etahat,n))
    a1 <- sum(ehat * mu(etahat + lam * ehat,n))
    a2 <- sum(ehat * mu(etahat - lam * ehat,n))
    a <- (a1 - 2 * a0 + a2)/(lam^2 * 6 * sighat^3)
    ## calculate bias bhat
    bvec <- numeric(p)
    eig <- eigen(S)
    evals <- (eig$values)^0.5
    evecs <- (eig$vectors)
    for(j in 1:p) {
        b1 <- tt(y + lambda * evals[j] * evecs[, j])
        b2 <- tt(y - lambda * evals[j] * evecs[, j])
        bvec[j] <- (b1 - 2 * thetahat + b2)/lambda^2
    }
    bhat <- sum(bvec)/2
    ## calculate quadratic coefficient cq
    delta <- dhat/sighat
    cq <- (tt(y + lambda * delta) - 2 * thetahat
           + tt(y - lambda * delta)) / (2 * sighat * lambda^2)
    ## calculate bias-correction constant z0
    curv <- bhat/sighat - cq
    z0 <- qnorm(2 * pnorm(a) * pnorm( - curv))	
    ## calculate Standard,ABC, and ABCq limits
    al <- c(alpha, rev(1 - alpha))
    za <- qnorm(al)
    z0a <- (z0 + za)/(1 - a * (z0 + za))
    z1a <- z0a + a * z0a^2	#   calculate endpoints
    standard <- thetahat + sighat * za
    ABC <- numeric(length(za))
    for(j in 1:length(za))
        ABC[j] <- tt(y + delta * z1a[j])
    ABCquad <- thetahat + sighat * (z1a + cq * z1a^2)
    limits <- cbind(al, ABC, ABCquad, standard)
    dimnames(limits) <- list(NULL, c("alpha", "ABC", "ABCquad",
                                     "Standard"))
    ## output in list form
    vl <- list(sys = syscall, 
              limits = limits, 
              stats =
                   list(thetahat=thetahat, sighat=sighat, bhat=bhat),
              constants = list(a=a, z0=z0, cq=cq),
              asym.05 = c(2 * a * 1.645, z0/1.645, cq * 1.645),
              call=call)
    vl$dhat <- dhat
    vl$ehat <- ehat
    return(vl)
}
