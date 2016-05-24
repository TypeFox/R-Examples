`nbind` <-
function (Yin, fm.X, t.i, fm.ga, gmat, nmat, totalit, rnull, betanull, gammanull, sigmakleinnull, psinull, Tau, alpha){
    d <- 0.005
    a <- 1
    asig <- 1
    bsig <- 0.005
    nmat <- as.matrix(nmat[, 2:length(nmat)])
    maxindex <- dim(nmat)[2]
    eigenv <- eigen(gmat - diag(nmat[, maxindex]))$values
    J <- dim(gmat)[1]
    Jasig <- J/2 + asig
    r <- rnull
    nb <- (dim(fm.X)[2] - 1)
    rvar <- 0.1
    psivar <- 0.1
    bvarvec <- cbind(0.05, 9, 7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
    gvarvec <- matrix(1, 1, J) * 1.4
    Yin <- as.matrix(Yin)
    fm.X <- as.matrix(fm.X)
    ny <- dim(Yin)[1]
    nx <- dim(fm.X)[2]
    fm.Xb <- as.matrix(fm.X[, 2:ncol(fm.X)])
    regindex <- fm.X[, 1]
    psi <- psinull
    invsigmasqnull <- 1/(sigmakleinnull^2)
    invsigmasq <- invsigmasqnull
    P0 <- diag(1/Tau^2, nb)
    bnull <- 0.1
    for (m in 1:1) {
        acceptg <- matrix(0, 1, length(gammanull))
        acceptb <- matrix(0, 1, length(betanull))
        acceptr <- double(1)
        acceptpsi <- double(1)
        betamat <- betanull
        gammamat <- gammanull
        psivec <- psinull
        invsigmasqnull <- 1/(sigmakleinnull^2)
        sigmavec <- invsigmasqnull
        rvec <- rnull
        bvec <- bnull
        be <- betanull
        ga <- gammanull
        r <- rnull
        b <- bnull
        psi <- psinull
        invsigmasq <- invsigmasqnull
        gammavec <- ga[regindex]
        mu <- t.i * exp(fm.Xb %*% be + gammavec)
        for (i in 0:100) {
            r_r <- .C("rnb", PACKAGE = "spatcounts", as.double(mu), 
                as.double(r), as.double(fm.X), as.integer(ny), 
                as.double(Yin), as.double(t.i), as.integer(acceptr), 
                as.double(rvar), as.double(a), as.double(b), 
                as.double(double(2)))[[11]]
            r <- r_r[1]
            acceptr <- r_r[2]
            rvec <- cbind(rvec, r)
            b <- rgamma(1, 1 + a, r + d)
            bvec <- cbind(bvec, b)
            beta_r <- .C("betanbind", PACKAGE = "spatcounts", 
                as.double(mu), as.integer(nb), as.double(be), 
                as.double(gammavec), as.double(r), as.double(fm.Xb), 
                as.integer(ny), as.double(Yin), as.double(t.i), 
                as.integer(acceptb), as.double(P0), as.double(double((nb + 
                  nb + ny))))[[12]]
            be <- beta_r[1:nb]
            acceptb <- beta_r[(nb + 1):(2 * nb)]
            mu <- beta_r[(2 * nb + 1):(2 * nb + ny)]
            betamat <- cbind(betamat, be)
            if (fm.ga == TRUE) {
                gamma_r <- .C("gammanbind", PACKAGE = "spatcounts", 
                  as.double(mu), as.integer(nb), as.double(be), 
                  as.integer(J), as.double(ga), as.double(r), 
                  as.double(psi), as.double(invsigmasq), as.double(fm.X), 
                  as.integer(ny), as.double(Yin), as.double(t.i), 
                  as.integer(acceptg), as.double(gvarvec), as.integer(nmat), 
                  as.integer(maxindex), as.double(double((J + 
                    J))))[[17]]
                ga <- gamma_r[1:J]
                acceptg <- gamma_r[(J + 1):(2 * J)]
                gammamat <- cbind(gammamat, ga)
                gammavec <- ga[regindex]
                mu <- t.i * exp(fm.Xb %*% be + gammavec)
                gKerng <- double(1)
                if (psi > 0) {
                  for (k in 1:J) {
                    gmatsum <- double(1)
                    for (j in 1:nmat[k, maxindex]) {
                      gmatsum <- gmatsum + ga[nmat[k, j]]
                    }
                    gKerng <- gKerng + ga[k] * (ga[k] * nmat[k, 
                      maxindex] - gmatsum)
                  }
                  gQg <- psi * gKerng + (t(ga) %*% ga)
                }
                else {
                  for (k in 1:J) {
                    gmatsum <- double(1)
                    for (j in 1:nmat[k, maxindex]) {
                      gmatsum <- gmatsum + ga[nmat[k, j]]
                    }
                    gKerng <- gKerng + ga[k] * (ga[k] * nmat[k, 
                      maxindex] + gmatsum)
                  }
                  gQg <- -psi * gKerng + (t(ga) %*% ga)
                }
                if (gQg > 0) {
                  invsigmasq <- rgamma(1, Jasig, 0.5 * gQg + 
                    bsig)
                  sigmavec <- cbind(sigmavec, invsigmasq)
                }
                psi_r <- .C("psimhbarbayern", PACKAGE = "spatcounts", 
                  as.double(psi), as.integer(J), as.double(ga), 
                  as.double(eigenv), as.double(invsigmasq), as.double(gKerng), 
                  as.integer(acceptpsi), as.double(psivar), as.double(alpha), 
                  as.double(double(2)))[[10]]
                psi <- psi_r[1]
                acceptpsi <- psi_r[2]
                psivec <- cbind(psivec, psi)
            }
            i <- i + 1
        }
        if (acceptr/102 > 0.6) {
            rvar <- rvar + rvar/6
        }
        if (acceptr/102 < 0.2) {
            rvar <- rvar - rvar/6
        }
        for (j in 1:length(be)) {
            if (acceptb[j]/102 > 0.6) {
                bvarvec[j] <- bvarvec[j] + bvarvec[j]/5
            }
            if (acceptb[j]/102 < 0.3) {
                bvarvec[j] <- bvarvec[j] - bvarvec[j]/5
            }
        }
        for (j in 1:length(ga)) {
            if (acceptg[j]/102 > 0.6) {
                gvarvec[j] <- gvarvec[j] + gvarvec[j]/10
            }
            if (acceptg[j]/102 < 0.3) {
                gvarvec[j] <- gvarvec[j] - gvarvec[j]/10
            }
        }
        if (acceptpsi/102 > 0.7) {
            psivar <- psivar + psivar/5
        }
        if (acceptpsi/102 < 0.3) {
            psivar <- psivar - psivar/5
        }
        m <- m + 1
    }
    betamat <- betanull
    gammamat <- gammanull
    psivec <- psinull
    invsigmasqnull <- 1/(sigmakleinnull^2)
    sigmavec <- invsigmasqnull
    rvec <- rnull
    bvec <- bnull
    invsigmasqnull <- 1/(sigmakleinnull^2)
    be <- betanull
    ga <- gammanull
    r <- rnull
    b <- bnull
    psi <- psinull
    invsigmasq <- invsigmasqnull
    acceptg <- matrix(0, 1, length(gammanull))
    acceptb <- matrix(0, 1, length(betanull))
    acceptr <- double(1)
    acceptpsi <- double(1)
    gammavec <- ga[regindex]
    mu <- t.i * exp(fm.Xb %*% be + gammavec)
    for (i in 0:totalit) {
        r_r <- .C("rnb", PACKAGE = "spatcounts", as.double(mu), 
            as.double(r), as.double(fm.X), as.integer(ny), as.double(Yin), 
            as.double(t.i), as.integer(acceptr), as.double(rvar), 
            as.double(a), as.double(b), as.double(double(2)))[[11]]
        r <- r_r[1]
        acceptr <- r_r[2]
        rvec <- cbind(rvec, r)
        b <- rgamma(1, 1 + a, r + d)
        bvec <- cbind(bvec, b)
        beta_r <- .C("betanbind", PACKAGE = "spatcounts", as.double(mu), 
            as.integer(nb), as.double(be), as.double(gammavec), 
            as.double(r), as.double(fm.Xb), as.integer(ny), as.double(Yin), 
            as.double(t.i), as.integer(acceptb), as.double(P0), 
            as.double(double((nb + nb + ny))))[[12]]
        be <- beta_r[1:nb]
        acceptb <- beta_r[(nb + 1):(2 * nb)]
        mu <- beta_r[(2 * nb + 1):(2 * nb + ny)]
        betamat <- cbind(betamat, be)
        if (fm.ga == TRUE) {
            gamma_r <- .C("gammanbind", PACKAGE = "spatcounts", 
                as.double(mu), as.integer(nb), as.double(be), 
                as.integer(J), as.double(ga), as.double(r), as.double(psi), 
                as.double(invsigmasq), as.double(fm.X), as.integer(ny), 
                as.double(Yin), as.double(t.i), as.integer(acceptg), 
                as.double(gvarvec), as.integer(nmat), as.integer(maxindex), 
                as.double(double((J + J))))[[17]]
            ga <- gamma_r[1:J]
            acceptg <- gamma_r[(J + 1):(2 * J)]
            gammamat <- cbind(gammamat, ga)
            gammavec <- ga[regindex]
            mu <- t.i * exp(fm.Xb %*% be + gammavec)
            gKerng <- double(1)
            if (psi > 0) {
                for (k in 1:J) {
                  gmatsum <- double(1)
                  for (j in 1:nmat[k, maxindex]) {
                    gmatsum <- gmatsum + ga[nmat[k, j]]
                  }
                  gKerng <- gKerng + ga[k] * (ga[k] * nmat[k, 
                    maxindex] - gmatsum)
                }
                gQg <- psi * gKerng + (t(ga) %*% ga)
            }
            else {
                for (k in 1:J) {
                  gmatsum <- double(1)
                  for (j in 1:nmat[k, maxindex]) {
                    gmatsum <- gmatsum + ga[nmat[k, j]]
                  }
                  gKerng <- gKerng + ga[k] * (ga[k] * nmat[k, 
                    maxindex] + gmatsum)
                }
                gQg <- -psi * gKerng + (t(ga) %*% ga)
            }
            if (gQg > 0) {
                invsigmasq <- rgamma(1, Jasig, 0.5 * gQg + bsig)
                sigmavec <- cbind(sigmavec, invsigmasq)
            }
            psi_r <- .C("psimhbarbayern", PACKAGE = "spatcounts", 
                as.double(psi), as.integer(J), as.double(ga), 
                as.double(eigenv), as.double(invsigmasq), as.double(gKerng), 
                as.integer(acceptpsi), as.double(psivar), as.double(alpha), 
                as.double(double(2)))[[10]]
            psi <- psi_r[1]
            acceptpsi <- psi_r[2]
            psivec <- cbind(psivec, psi)
        }
        i <- i + 1
    }
    minag <- min(acceptg)
    maxag <- max(acceptg)
    cat("acceptb/(i+1)  ", acceptb/(i + 1))
    cat("\n")
    cat("acceptga1/i acceptga2/(i+1)  ", cbind(minag/(i + 1), 
        maxag/(i + 1)))
    cat("\n")
    cat("acceptr/(i+1) ", acceptr/(i + 1))
    cat("\n")
    cat("acceptpsi/(i+1)  ", acceptpsi/(i + 1))
    cat("\n")
    range.gamma <- c(minag/(i + 1), maxag/(i + 1))
    if (fm.ga == FALSE) {
        Coefficients <- length(betanull) + length(rnull)
        gammamat <- matrix(0, J, i + 1)
        sigmavec <- matrix(0, 1, i + 1)
        psivec <- matrix(0, 1, i + 1)
    }
    else {
        Coefficients <- length(betanull) + length(rnull) + length(gmat)
    }
    nb.data <- list(acceptb = acceptb/(i + 1), acceptga = range.gamma, 
        acceptpsi = acceptpsi/(i + 1), acceptr = acceptr/(i + 
            1), beta = betamat, gamma = gammamat, invsigsq = sigmavec, 
        psi = psivec, r = rvec, b = bvec, Coefficients = Coefficients, 
        t.i = t.i)
    return(nb.data)
}

