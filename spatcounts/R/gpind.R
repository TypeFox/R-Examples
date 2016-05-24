`gpind` <-
function (Yin, fm.X, t.i, fm.ga = TRUE, gmat, nmat, totalit, phinull, betanull, gammanull, sigmakleinnull, psinull, Tau, alpha){
    asig <- 1
    bsig <- 0.005
    nmat <- as.matrix(nmat[, 2:length(nmat)])
    maxindex <- ncol(nmat)
    eigenv <- eigen(gmat - diag(nmat[, maxindex]))$values
    J <- ncol(gmat)
    Jasig <- J/2 + asig
    phi <- phinull
    nb <- (dim(fm.X)[2] - 1)
    phivar <- 1.05
    psivar <- 0.1
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
    for (m in 1:1) {
        acceptg <- matrix(0, 1, length(gammanull))
        acceptb <- matrix(0, 1, length(betanull))
        acceptphi <- double(1)
        acceptpsi <- double(1)
        phivec <- phi
        betamat <- betanull
        gammamat <- gammanull
        psivec <- psinull
        sigmavec <- invsigmasqnull
        be <- betanull
        ga <- gammanull
        phi <- phinull
        psi <- psinull
        invsigmasq <- invsigmasqnull
        gammavec <- ga[regindex]
        mu <- t.i * exp(fm.Xb %*% be + gammavec)
        for (i in 0:100) {
            phi_r <- .C("phigpind", PACKAGE = "spatcounts", as.double(phi), 
                as.double(mu), as.double(Yin), as.integer(acceptphi), 
                as.integer(ny), as.double(double(2)))[[6]]
            phi <- phi_r[1]
            acceptphi <- phi_r[2]
            phivec <- cbind(phivec, phi)
            beta_r <- .C("betaindbisection", PACKAGE = "spatcounts", 
                as.double(mu), as.integer(nb), as.double(be), 
                as.double(phi), as.integer(ny), as.double(gammavec), 
                as.double(fm.Xb), as.double(Yin), as.double(t.i), 
                as.double(P0), as.integer(acceptb), as.double(double((nb + 
                  nb + ny))))[[12]]
            be <- beta_r[1:nb]
            acceptb <- beta_r[(nb + 1):(2 * nb)]
            mu <- beta_r[(2 * nb + 1):(2 * nb + ny)]
            betamat <- cbind(betamat, be)
            if (fm.ga == TRUE) {
                gamma_r <- .C("gammaindbisection", PACKAGE = "spatcounts", 
                  as.integer(nb), as.double(be), as.integer(J), 
                  as.double(ga), as.integer(ny), as.double(mu), 
                  as.double(phi), as.double(psi), as.double(invsigmasq), 
                  as.double(fm.X), as.double(Yin), as.double(t.i), 
                  as.integer(acceptg), as.integer(nmat), as.integer(maxindex), 
                  as.double(double((J + J))))[[16]]
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
        if (acceptphi/102 > 0.6) {
            phivar <- phivar + phivar/5
        }
        if (acceptphi/102 < 0.3) {
            phivar <- phivar - phivar/5
        }
        if (acceptpsi/102 > 0.6) {
            psivar <- psivar + psivar/5
        }
        if (acceptpsi/102 < 0.3) {
            psivar <- psivar - psivar/5
        }
        m <- m + 1
    }
    phivec <- phinull
    betamat <- betanull
    gammamat <- gammanull
    psivec <- psinull
    sigmavec <- invsigmasqnull
    be <- betanull
    ga <- gammanull
    phi <- phinull
    psi <- psinull
    invsigmasq <- invsigmasqnull
    acceptg <- matrix(0, 1, length(gammanull))
    acceptb <- matrix(0, 1, length(betanull))
    acceptphi <- double(1)
    acceptpsi <- double(1)
    gammavec <- ga[regindex]
    mu <- t.i * exp(fm.Xb %*% be + gammavec)
    for (i in 0:totalit) {
        phi_r <- .C("phigpind", PACKAGE = "spatcounts", as.double(phi), 
            as.double(mu), as.double(Yin), as.integer(acceptphi), 
            as.integer(ny), as.double(double(2)))[[6]]
        phi <- phi_r[1]
        acceptphi <- phi_r[2]
        phivec <- cbind(phivec, phi)
        beta_r <- .C("betaindbisection", PACKAGE = "spatcounts", 
            as.double(mu), as.integer(nb), as.double(be), as.double(phi), 
            as.integer(ny), as.double(gammavec), as.double(fm.Xb), 
            as.double(Yin), as.double(t.i), as.double(P0), as.integer(acceptb), 
            as.double(double((nb + nb + ny))))[[12]]
        be <- beta_r[1:nb]
        acceptb <- beta_r[(nb + 1):(2 * nb)]
        mu <- beta_r[(2 * nb + 1):(2 * nb + ny)]
        betamat <- cbind(betamat, be)
        if (fm.ga == TRUE) {
            gamma_r <- .C("gammaindbisection", PACKAGE = "spatcounts", 
                as.integer(nb), as.double(be), as.integer(J), 
                as.double(ga), as.integer(ny), as.double(mu), 
                as.double(phi), as.double(psi), as.double(invsigmasq), 
                as.double(fm.X), as.double(Yin), as.double(t.i), 
                as.integer(acceptg), as.integer(nmat), as.integer(maxindex), 
                as.double(double((J + J))))[[16]]
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
    cat("acceptb/(i+1) ", acceptb/(i + 1))
    cat("\n")
    cat("acceptga1/i acceptga2/(i+1)  ", cbind(minag/(i + 1), 
        maxag/(i + 1)))
    cat("\n")
    cat("acceptphi/(i+1) ", acceptphi/(i + 1))
    cat("\n")
    cat("acceptpsi/(i+1)  ", acceptpsi/(i + 1))
    cat("\n")
    range.gamma <- c(minag/(i + 1), maxag/(i + 1))
    if (fm.ga == FALSE) {
        Coefficients <- length(betanull) + length(phinull)
        gammamat <- matrix(0, J, i + 1)
        sigmavec <- matrix(0, 1, i + 1)
        psivec <- matrix(0, 1, i + 1)
    }
    else {
        Coefficients <- length(betanull) + length(phinull) + 
            length(gmat)
    }
    gp.data <- list(acceptb = acceptb/(i + 1), acceptga = range.gamma, 
        acceptpsi = acceptpsi/(i + 1), acceptphi = acceptphi/(i + 
            1), beta = betamat, gamma = gammamat, invsigsq = sigmavec, 
        psi = psivec, phi = phivec, Coefficients = Coefficients, 
        t.i = t.i)
    return(gp.data)
}

