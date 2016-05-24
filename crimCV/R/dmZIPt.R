dmZIPt <-
function (Dat, X, ng, rcv = TRUE, init = 20, Risk = NULL) 
{
    initpr <- function(Dat, X, Risk, ni, no, npp, ng, npop) {
        pparam <- rep(0, (npp + 2) * ng * npop)
        pllike <- rep(0, npop)
        Datm <- Dat
        Datm[Datm < 0] <- 0
        Dmu <- apply(Datm, 2, mean)
        Dmu <- rep(Dmu, ni)
        Dmu <- matrix(Dmu, no, ni)
        Dmu <- t(Dmu)
        Dmu[Dat > -0.1] <- 0
        Datm <- Datm + Dmu
        Frtr <- .Fortran("r_dmzipt_init_param", as.double(X), 
            as.double(Datm), as.double(Risk), pparam = as.double(pparam), 
            pllike = as.double(pllike), as.integer(ni), as.integer(no), 
            as.integer(npp), as.integer(ng), as.integer(npop))
        out <- NULL
        out$pparam <- matrix(Frtr$pparam, npop, (npp + 2) * ng)
        out$pllike <- Frtr$pllike
        out
    }
    ni <- nrow(Dat)
    no <- ncol(Dat)
    npp <- ncol(X)
    if (nrow(X) != no) 
        stop("dmZIPt: Incorrect input nrows(X) must equal ncol(Dat)!")
    if (is.null(Risk)) {
        Risk <- matrix(1, ni, no)
    }
    nn <- ni * no
    npr <- (npp + 1) * ng + ng
    cat("|----Initialize----:")
    flush.console()
    npop <- init * ng
    outi <- initpr(Dat, X, Risk, ni, no, npp, ng, npop)
    if (is.null(outi)) {
        cat(" [FAILED] \n")
        return(outi)
    }
    cat(" [DONE] \n")
    flush.console()
    pparam <- matrix(outi$pparam, npop, (npp + 2) * ng)
    bparam <- NULL
    bllike <- NULL
    llike <- -1e+20
    cat("|-----Fitting------:")
    flush.console()
    for (i in 1:npop) {
        param <- pparam[i, ]
        param <- matrix(param, npp + 2, ng)
        beta <- param[1:npp, ]
        tau <- param[npp + 1, ]
        prob <- param[npp + 2, ]
        ggt <- matrix(0, ni, ng)
        Info <- matrix(0, npp * ng + 2 * ng - 1, npp * ng + 2 * 
            ng - 1)
        tFrtr <- .Fortran("r_dmzipt", as.double(X), as.double(Dat), 
            as.double(Risk), ggt = as.double(ggt), prob = as.double(prob), 
            beta = as.double(beta), tau = as.double(tau), llike = as.double(0), 
            Info = as.double(Info), as.integer(nn), as.integer(ni), 
            as.integer(no), as.integer(npp), as.integer(ng), 
            err = as.integer(0))
        if (tFrtr$err != 0) 
            next
        if (tFrtr$llike > llike) {
            Frtr <- tFrtr
            llike <- tFrtr$llike
        }
        if (i == 1) {
            beta <- matrix(tFrtr$beta, npp, ng)
            tau <- exp(tFrtr$tau)
            prob <- tFrtr$prob
            param <- rbind(beta, tau, prob)
            param <- c(param)
            bparam <- cbind(bparam, param)
            bllike <- c(bllike, tFrtr$llike)
        }
        else {
            if (!(min(abs(tFrtr$llike - bllike)) < 1e-05)) {
                beta <- matrix(tFrtr$beta, npp, ng)
                tau <- exp(tFrtr$tau)
                prob <- tFrtr$prob
                param <- rbind(beta, tau, prob)
                param <- c(param)
                bparam <- cbind(bparam, param)
                bllike <- c(bllike, tFrtr$llike)
            }
        }
    }
    if (length(bllike) > 1) {
        rllike <- order(-bllike)
        bllike <- bllike[rllike]
        bparam <- bparam[, rllike]
    }
    out <- NULL
    if (Frtr$err != 0) {
        cat(" [FAILED] \n")
        return(out)
    }
    cat(" [DONE] \n")
    flush.console()
    out$beta <- matrix(Frtr$beta, npp, ng)
    out$tau <- exp(Frtr$tau)
    out$prob <- Frtr$prob
    out$gwt <- matrix(Frtr$ggt, ni, ng)
    out$Info <- matrix(Frtr$Info, npp * ng + 2 * ng - 1, npp * 
        ng + 2 * ng - 1)
    llike <- Frtr$llike
    out$llike <- llike
    out$AIC <- -2 * llike + 2 * (ng * (npp + 1) + ng - 1)
    out$BIC <- -2 * llike + (ng * (npp + 1) + ng - 1) * log(nn)
    out$bparam <- bparam
    out$bllike <- bllike
    out$ni <- ni
    out$ng <- ng
    out$X <- X
    if (rcv) {
        cat("|---CV/Jackknife---: \n")
        flush.console()
        param <- matrix(0, npr, ni)
        parm <- c(c(out$beta), out$tau, out$prob)
        iseq <- 1:ni
        imask <- rep(TRUE, ni)
        nn <- (ni - 1) * no
        cvllike <- rep(0, ni)
        cvDat <- matrix(0, ni, no)
        cvDat2 <- matrix(0, ni, no)
        Info <- matrix(0, npp * ng + 2 * ng - 1, npp * ng + 2 * 
            ng - 1)
        errcnt <- 0
        for (i in 1:ni) {
            cat("                   : i =", i, "\n")
            flush.console()
            imask[i] <- FALSE
            indi <- iseq[imask]
            tDat <- Dat[indi, ]
            tRisk <- Risk[indi, ]
            ggt <- out$gwt[indi, ]
            beta <- out$beta
            tau <- log(out$tau)
            if (ng > 1) {
                prob <- out$prob
            }
            else {
                prob <- 1
            }
            cFrtr <- .Fortran("r_dmzipt", as.double(X), as.double(tDat), 
                as.double(tRisk), ggt = as.double(ggt), prob = as.double(prob), 
                beta = as.double(beta), tau = as.double(tau), 
                llike = as.double(0), Info = as.double(Info), 
                as.integer(nn), as.integer(ni - 1), as.integer(no), 
                as.integer(npp), as.integer(ng), err = as.integer(0))
            if (cFrtr$err != 0) {
                param[, i] <- c(out$beta, out$tau, out$prob)
                cvDat[i, ] <- NA
                errcnt <- errcnt + 1
                imask[i] <- TRUE
                next
            }
            beta <- cFrtr$beta
            tau <- exp(cFrtr$tau)
            prob <- cFrtr$prob
            param[, i] <- c(beta, tau, prob)
            beta <- matrix(beta, npp, ng)
            Xb <- X %*% beta
            Zg <- -tau * t(Xb)
            Xb <- log(Risk[i, ]) + Xb
            Zg <- t(Zg)
            lam <- exp(Xb)
            qq <- exp(Zg)
            qq <- qq/(1 + qq)
            mu <- (1 - qq) * lam
            mu2 <- mu
            llike = log(exp(Zg) + exp(-exp(Xb)))
            DD <- Dat[i, ]
            DD[DD < 0] <- 0
            llike[DD > 0.5, ] <- 0
            tllike = c(DD) * Xb - exp(Xb) - lgamma(DD + 1)
            tllike[DD < 0.5, ] <- 0
            llike = llike + tllike
            llike = llike - log(1 + exp(Zg))
            llike[Dat[i, ] < 0, ] <- 0
            ggt <- apply(llike, 2, sum)
            ggt <- log(prob) + ggt
            rscl <- max(ggt)
            ggt <- exp(ggt - rscl)
            cvllike[i] <- sum(ggt)
            ggt <- ggt/cvllike[i]
            cvllike[i] <- log(cvllike[i]) + rscl
            for (j in 1:ng) {
                mu[, j] <- ggt[j] * mu[, j]
                mu2[, j] <- prob[j] * mu[, j]
            }
            cvDat[i, ] <- apply(mu, 1, sum)
            cvDat2[i, ] <- apply(mu2, 1, sum)
            imask[i] <- TRUE
        }
        cat("|---CV/Jackknife---: [DONE] \n")
        flush.console()
        cv <- mean(abs(Dat - cvDat), na.rm = TRUE)
        cvunc <- mean(abs(Dat - cvDat2), na.rm = TRUE)
        cv2 <- mean((Dat - cvDat)^2, na.rm = TRUE)
        dparm <- ni * parm - (ni - 1) * param
        jVar <- var(t(dparm))/ni
        out$cv <- cv
        out$cvunc <- cvunc
        out$cv2 <- cv2
        out$lcv <- -2 * sum(cvllike)
        out$param <- param
        out$jVar <- jVar
        out$errcnt <- errcnt
    }
    class(out) <- c("dmZIPt", class(out))
    out
}
