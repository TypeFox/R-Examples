`zipind` <-
function (Yin, fm.X, t.i, fm.ga, gmat, nmat, totalit, omeganull, betanull, gammanull, sigmakleinnull, psinull, Tau, alpha){
    asig <- 1
    bsig <- 0.005
    nmat <- as.matrix(nmat[, 2:length(nmat)])
    maxindex <- ncol(nmat)
    eigenv <- eigen(gmat - diag(nmat[, maxindex]))$values
    J <- ncol(gmat)
    Jasig <- J/2 + asig
    phinull <- 1
    nb <- (dim(fm.X)[2] - 1)
    psivar <- 0.1
    omegavar <- 0.07
    bvarvec <- c(0.1, 8, 7, 0.1, 0.07, 0.2, 0.1, 0.1, 0.08)
    gvarvec <- matrix(1, 1, J) * 0.4
    nyzero <- length(Yin[Yin == 0])
    znull <- rbinom(nyzero, 1, omeganull)
    Yin <- as.matrix(Yin)
    fm.X <- as.matrix(fm.X)
    ny <- dim(Yin)[1]
    yn1 <- length(Yin[Yin > 0])
    nx <- dim(fm.X)[2]
    fm.Xb <- as.matrix(fm.X[, 2:ncol(fm.X)])
    regindex <- fm.X[, 1]
    psi <- psinull
    phi <- phinull
    invsigmasqnull <- 1/(sigmakleinnull^2)
    invsigmasq <- invsigmasqnull
    P0 <- diag(1/Tau^2, nb)
    zfull <- matrix(0, 1, ny)
    zfull[Yin == 0] <- znull
    data1 <- cbind(Yin, fm.X, t.i)
    z0data <- data1[zfull == 0, ]
    for (m in 1:1) {
        acceptb <- matrix(0, 1, length(betanull))
        acceptg <- matrix(0, 1, length(gammanull))
        acceptpsi <- double(1)
        acceptomega <- double(1)
        acceptbint <- double(1)
        gammamat <- gammanull
        psivec <- psinull
        sigmavec <- invsigmasqnull
        ga <- gammanull
        betamat <- betanull
        zvec <- znull
        be <- betanull
        z <- znull
        omega <- omeganull
        omegavec <- omega
        gammavec <- ga[regindex]
        mu <- t.i * exp(fm.Xb %*% be + gammavec)
        for (i in 0:100) {
            beta_r <- .C("betaintercept", PACKAGE = "spatcounts", 
                as.double(mu), as.integer(nb), as.double(be), 
                as.double(phi), as.integer(ny), as.double(gammavec), 
                as.double(omega), as.double(fm.Xb), as.double(Yin), 
                as.double(t.i), as.double(P0), as.integer(acceptbint), 
                as.double(double((nb + nb + ny))))[[13]]
            be <- beta_r[1:nb]
            acceptbint <- beta_r[nb + 1]
            mu <- beta_r[(2 * nb + 1):(2 * nb + ny)]
            omega_r <- .C("omegazigpind", PACKAGE = "spatcounts", 
                as.double(omega), as.double(mu), as.double(phi), 
                as.double(fm.X), as.double(Yin), as.double(t.i), 
                as.integer(acceptomega), as.integer(ny), as.integer(yn1), 
                as.double(double(2)))[[10]]
            omega <- omega_r[1]
            acceptomega <- omega_r[2]
            omegavec <- cbind(omegavec, omega)
            mu0 <- mu[Yin == 0]
            z <- rbinom(nyzero, 1, omega/(omega + (1 - omega) * 
                exp(-mu0)))
            zvec <- zvec + z
            zfull <- matrix(0, 1, ny)
            zfull[Yin == 0] <- z
            beta_r <- .C("betazigpindbisection", PACKAGE = "spatcounts", 
                as.double(mu), as.integer(nb), as.double(be), 
                as.double(phi), as.integer(ny), as.double(gammavec), 
                as.double(fm.X), as.integer(zfull), as.double(Yin), 
                as.double(t.i), as.double(P0), as.integer(acceptb), 
                as.double(double((nb + nb + ny))))[[13]]
            be <- beta_r[1:nb]
            acceptb <- beta_r[(nb + 1):(2 * nb)]
            mu <- beta_r[(2 * nb + 1):(2 * nb + ny)]
            betamat <- cbind(betamat, be)
            if (fm.ga == TRUE) {
                gamma_r <- .C("gammazigpindbisection", PACKAGE = "spatcounts", 
                  as.integer(nb), as.double(be), as.integer(J), 
                  as.double(ga), as.integer(ny), as.double(mu), 
                  as.double(phi), as.double(psi), as.double(invsigmasq), 
                  as.double(fm.X), as.double(Yin), as.integer(zfull), 
                  as.double(t.i), as.integer(acceptg), as.integer(nmat), 
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
        if (acceptpsi/102 > 0.7) {
            psivar <- psivar + psivar/5
        }
        if (acceptpsi/102 < 0.3) {
            psivar <- psivar - psivar/5
        }
        for (j in 1:length(ga)) {
            if (acceptg[j]/102 > 0.6) {
                gvarvec[j] <- gvarvec[j] + gvarvec[j]/10
            }
            if (acceptg[j]/102 < 0.3) {
                gvarvec[j] <- gvarvec[j] - gvarvec[j]/10
            }
        }
        m <- m + 1
    }
    acceptg <- matrix(0, 1, length(gammanull))
    acceptb <- matrix(0, 1, length(betanull))
    acceptbint <- double(1)
    acceptpsi <- double(1)
    acceptomega <- double(1)
    gammamat <- gammanull
    betamat <- betanull
    zvec <- znull
    psivec <- psinull
    invsigmasqnull <- 1/(sigmakleinnull^2)
    sigmavec <- invsigmasqnull
    ga <- gammanull
    psi <- psinull
    invsigmasq <- invsigmasqnull
    be <- betanull
    z <- znull
    omega <- omeganull
    omegavec <- omega
    gammavec <- ga[regindex]
    mu <- t.i * exp(fm.Xb %*% be + gammavec)
    for (i in 0:totalit) {
        beta_r <- .C("betaintercept", PACKAGE = "spatcounts", 
            as.double(mu), as.integer(nb), as.double(be), as.double(phi), 
            as.integer(ny), as.double(gammavec), as.double(omega), 
            as.double(fm.Xb), as.double(Yin), as.double(t.i), 
            as.double(P0), as.integer(acceptbint), as.double(double((nb + 
                nb + ny))))[[13]]
        be <- beta_r[1:nb]
        acceptbint <- beta_r[nb + 1]
        mu <- beta_r[(2 * nb + 1):(2 * nb + ny)]
        omega_r <- .C("omegazigpind", PACKAGE = "spatcounts", 
            as.double(omega), as.double(mu), as.double(phi), 
            as.double(fm.X), as.double(Yin), as.double(t.i), 
            as.integer(acceptomega), as.integer(ny), as.integer(yn1), 
            as.double(double(2)))[[10]]
        omega <- omega_r[1]
        acceptomega <- omega_r[2]
        omegavec <- cbind(omegavec, omega)
        mu0 <- mu[Yin == 0]
        z <- rbinom(nyzero, 1, omega/(omega + (1 - omega) * exp(-mu0)))
        zvec <- zvec + z
        zfull <- matrix(0, 1, ny)
        zfull[Yin == 0] <- z
        beta_r <- .C("betazigpindbisection", PACKAGE = "spatcounts", 
            as.double(mu), as.integer(nb), as.double(be), as.double(phi), 
            as.integer(ny), as.double(gammavec), as.double(fm.X), 
            as.integer(zfull), as.double(Yin), as.double(t.i), 
            as.double(P0), as.integer(acceptb), as.double(double((nb + 
                nb + ny))))[[13]]
        be <- beta_r[1:nb]
        acceptb <- beta_r[(nb + 1):(2 * nb)]
        mu <- beta_r[(2 * nb + 1):(2 * nb + ny)]
        betamat <- cbind(betamat, be)
        if (fm.ga == TRUE) {
            gamma_r <- .C("gammazigpindbisection", PACKAGE = "spatcounts", 
                as.integer(nb), as.double(be), as.integer(J), 
                as.double(ga), as.integer(ny), as.double(mu), 
                as.double(phi), as.double(psi), as.double(invsigmasq), 
                as.double(fm.X), as.double(Yin), as.integer(zfull), 
                as.double(t.i), as.integer(acceptg), as.integer(nmat), 
                as.integer(maxindex), as.double(double((J + J))))[[17]]
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
    acceptbeta <- c(acceptbint, acceptb[2:nb])
    cat("acceptb/(i+1)  ", acceptbeta/(i+1))
    cat("\n")
    cat("acceptga1/i acceptga2/(i+1)  ", cbind(minag/(i + 1), 
        maxag/(i + 1)))
    cat("\n")
    cat("acceptomega/(i+1) ", acceptomega/(i + 1))
    cat("\n")
    cat("acceptpsi/(i+1)  ", acceptpsi/(i + 1))
    cat("\n")
    range.gamma <- c(minag/(i + 1), maxag/(i + 1))
    phivec <- matrix(1, 1, i + 1)
    if (fm.ga == FALSE) {
        Coefficients <- length(betanull) + length(omeganull)
        gammamat <- matrix(0, J, i + 1)
        sigmavec <- matrix(0, 1, i + 1)
        psivec <- matrix(0, 1, i + 1)
    }
    else {
        Coefficients <- length(betanull) + length(omeganull) + 
            length(gmat)
    }
    zip.data <- list(acceptb = acceptbeta/(i+1), acceptga = range.gamma, 
        acceptpsi = acceptpsi/(i + 1), acceptomega = acceptomega/(i + 
            1), beta = betamat, gamma = gammamat, invsigsq = sigmavec, 
        psi = psivec, phi = phivec, omega = omegavec, z = zvec, 
        Coefficients = Coefficients, t.i = t.i)
    return(zip.data)
}

