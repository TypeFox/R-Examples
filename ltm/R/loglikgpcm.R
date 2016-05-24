loglikgpcm <-
function (thetas, constraint) {
    betas <- betas.gpcm(thetas, p, ncatg, constraint)
    log.crf <- crf.GPCM(betas, Z, IRT.param, log = TRUE)
    log.p.xz <- matrix(0, nfreqs, length(Z))
    for (j in 1:p) {
        log.pr <- log.crf[[j]]
        xj <- X[, j]
        na.ind <- is.na(xj)
        log.pr <- log.pr[xj, ]
        if (any(na.ind))
            log.pr[na.ind, ] <- 0
        log.p.xz <- log.p.xz + log.pr
    }
    p.x <- rep(exp(log.p.xz) %*% GHw, obs)
    - sum(log(p.x))
}
