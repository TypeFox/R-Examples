dmZIP <-
function (Dat, X, Z, ng, rcv = TRUE, init = 20, Risk = NULL) 
{
    initpr <- function(Dat, X, Z, Risk, ni, no, npp, npl, ng, 
        npop) {
        pparam <- rep(0, (npp + npl + 1) * ng * npop)
        pllike <- rep(0, npop)
        Datm <- Dat
        Datm[Datm < 0] <- 0
        Dmu <- apply(Datm, 2, mean)
        Dmu <- rep(Dmu, ni)
        Dmu <- matrix(Dmu, no, ni)
        Dmu <- t(Dmu)
        Dmu[Dat > -0.1] <- 0
        Datm <- Datm + Dmu
        Frtr <- .Fortran("r_dmzip_init_param", as.double(X), 
            as.double(Z), as.double(Datm), as.double(Risk), pparam = as.double(pparam), 
            pllike = as.double(pllike), as.integer(ni), as.integer(no), 
            as.integer(npp), as.integer(npl), as.integer(ng), 
            as.integer(npop))
        out <- NULL
        out$pparam <- matrix(Frtr$pparam, npop, (npp + npl + 
            1) * ng)
        out$pllike <- Frtr$pllike
        out
    }
    ni <- nrow(Dat)
    no <- ncol(Dat)
    npp <- ncol(X)
    npl <- ncol(Z)
    if (nrow(X) != nrow(Z)) 
        stop("dmZIP: Incorrect input nrows(X) must equal ncol(Z)!")
    if (nrow(X) != no) 
        stop("dmZIP: Incorrect input nrows(X) must equal ncol(Dat)!")
    if (is.null(Risk)) {
        Risk <- matrix(1, ni, no)
    }
    nn <- ni * no
    npr <- (npp + npl) * ng + ng
    cat("|----Initialize----:")
    flush.console()
    npop <- init * ng
    outi <- initpr(Dat, X, Z, Risk, ni, no, npp, npl, ng, npop)
    if (is.null(outi)) {
        cat(" [FAILED] \n")
        return(outi)
    }
    cat(" [DONE] \n")
    flush.console()
    pparam <- matrix(outi$pparam, npop, (npp + npl + 1) * ng)
    bparam <- NULL
    bllike <- NULL
    llike <- -1e+20
    cat("|-----Fitting------:")
    flush.console()
    for (i in 1:npop) {
        param <- pparam[i, ]
        param <- matrix(param, npp + npl + 1, ng)
        beta <- param[1:npp, ]
        gamma <- param[(npp + 1):(npp + npl), ]
        prob <- param[npp + npl + 1, ]
        ggt <- matrix(0, ni, ng)
        Info <- matrix(0, (npp + npl) * ng + ng - 1, (npp + npl) * 
            ng + ng - 1)
        tFrtr <- .Fortran("r_dmzip", as.double(X), as.double(Z), 
            as.double(Dat), as.double(Risk), ggt = as.double(ggt), 
            prob = as.double(prob), beta = as.double(beta), gamma = as.double(gamma), 
            llike = as.double(0), Info = as.double(Info), as.integer(nn), 
            as.integer(ni), as.integer(no), as.integer(npp), 
            as.integer(npl), as.integer(ng), err = as.integer(0))
        if (tFrtr$err != 0) 
            next
        if (tFrtr$llike > llike) {
            Frtr <- tFrtr
            llike <- tFrtr$llike
        }
        if (i == 1) {
            beta <- matrix(tFrtr$beta, npp, ng)
            gamma <- matrix(tFrtr$gamma, npl, ng)
            prob <- tFrtr$prob
            param <- rbind(beta, gamma, prob)
            param <- c(param)
            bparam <- cbind(bparam, param)
            bllike <- c(bllike, tFrtr$llike)
        }
        else {
            if (!(min(abs(tFrtr$llike - bllike)) < 1e-05)) {
                beta <- matrix(tFrtr$beta, npp, ng)
                gamma <- matrix(tFrtr$gamma, npl, ng)
                prob <- tFrtr$prob
                param <- rbind(beta, gamma, prob)
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
    out$gamma <- matrix(Frtr$gamma, npl, ng)
    out$prob <- Frtr$prob
    out$gwt <- matrix(Frtr$ggt, ni, ng)
    out$Info <- matrix(Frtr$Info, (npp + npl) * ng + ng - 1, 
        (npp + npl) * ng + ng - 1)
    llike <- Frtr$llike
    out$llike <- llike
    out$AIC <- -2 * llike + 2 * (ng * (npp + npl) + ng - 1)
    out$BIC <- -2 * llike + (ng * (npp + npl) + ng - 1) * log(nn)
    out$bparam <- bparam
    out$bllike <- bllike
    out$ni <- ni
    out$ng <- ng
    out$X <- X
    out$Z <- Z
    if (rcv) {
        cat("|---CV/Jackknife---: \n")
        flush.console()
        param <- matrix(0, npr, ni)
        parm <- c(c(out$beta), c(out$gamma), out$prob)
        iseq <- 1:ni
        imask <- rep(TRUE, ni)
        nn <- (ni - 1) * no
        cvllike <- rep(0, ni)
        cvDat <- matrix(0, ni, no)
        Info <- matrix(0, (npp + npl) * ng + ng - 1, (npp + npl) * 
            ng + ng - 1)
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
            gamma <- out$gamma
            if (ng > 1) {
                prob <- out$prob
            }
            else {
                prob <- 1
            }
            cFrtr <- .Fortran("r_dmzip", as.double(X), as.double(Z), 
                as.double(tDat), as.double(tRisk), ggt = as.double(ggt), 
                prob = as.double(prob), beta = as.double(beta), 
                gamma = as.double(gamma), llike = as.double(0), 
                Info = as.double(Info), as.integer(nn), as.integer(ni - 
                  1), as.integer(no), as.integer(npp), as.integer(npl), 
                as.integer(ng), err = as.integer(0))
            if (cFrtr$err != 0) {
                param[, i] <- c(out$beta, out$gamma, out$prob)
                cvDat[i, ] <- NA
                errcnt <- errcnt + 1
                imask[i] <- TRUE
                next
            }
            beta <- cFrtr$beta
            gamma <- cFrtr$gamma
            prob <- cFrtr$prob
            param[, i] <- c(beta, gamma, prob)
            beta <- matrix(beta, npp, ng)
            gamma <- matrix(gamma, npl, ng)
            Xb <- X %*% beta
            Zg <- Z %*% gamma
            Xb <- log(Risk[i, ]) + Xb
            lam <- exp(Xb)
            qq <- exp(Zg)
            qq <- qq/(1 + qq)
            mu <- (1 - qq) * lam
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
            }
            cvDat[i, ] <- apply(mu, 1, sum)
            imask[i] <- TRUE
        }
        cat("|---CV/Jackknife---: [DONE] \n")
        flush.console()
        cv <- mean(abs(Dat - cvDat), na.rm = TRUE)
        cv2 <- mean((Dat - cvDat)^2, na.rm = TRUE)
        dparm <- ni * parm - (ni - 1) * param
        jVar <- var(t(dparm))/ni
        out$cv <- cv
        out$cv2 <- cv2
        out$lcv <- -2 * sum(cvllike)
        out$param <- param
        out$jVar <- jVar
        out$errcnt <- errcnt
    }
    class(out) <- c("dmZIP", class(out))
    out
}
